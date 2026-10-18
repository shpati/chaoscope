/* chaoscope_gui.c -- Win32 front end for chaoscope.exe
 *
 * Every picture is produced by running chaoscope.exe as a child process; the
 * BMP it writes is loaded and shown.  Systems, their parameter names and
 * defaults, and the palettes are read from `chaoscope.exe --list-tsv` at
 * start-up, so systems added to chaoscope.c show up here without touching
 * this file.
 *
 *   change any control   quick low-resolution preview, then a full-quality
 *                        render once nothing has changed for ~1 s
 *   drag on the image    rotate (roll / yaw)      right-drag: pitch
 *   mouse wheel          zoom                      Spin: continuous yaw
 *
 * Build (MinGW-w64):
 *   windres chaoscope_gui.rc -O coff -o chaoscope_gui_res.o
 *   gcc -O2 -mwindows -o chaoscope_gui.exe chaoscope_gui.c chaoscope_gui_res.o \
 *       -lcomctl32 -lcomdlg32 -lgdi32 -luser32
 * Build (MSVC, Developer Command Prompt):
 *   rc chaoscope_gui.rc
 *   cl /O2 chaoscope_gui.c chaoscope_gui.res /link /subsystem:windows
 * Usage:  chaoscope_gui.exe [path\to\chaoscope.exe]
 *         (default: chaoscope.exe next to this program, then PATH)
 */

#ifndef WINVER
#define WINVER 0x0501
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0501
#endif
#ifndef _WIN32_IE
#define _WIN32_IE 0x0600
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <commctrl.h>
#include <commdlg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#ifdef _MSC_VER
#pragma comment(lib, "user32.lib")
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "comctl32.lib")
#pragma comment(lib, "comdlg32.lib")
#if _MSC_VER < 1900
#define snprintf _snprintf
#define vsnprintf _vsnprintf
#endif
#endif

#ifndef ICC_STANDARD_CLASSES
#define ICC_STANDARD_CLASSES 0x00004000
#endif

/* ------------------------------------------------------------------ */
/* Constants                                                           */
/* ------------------------------------------------------------------ */

#define MAXSYS          128
#define MAXPAL          16
#define UIP             6        /* parameter fields offered per system   */
#define PREVIEW_MAX_DIM 800      /* longest preview side in pixels        */
#define PREVIEW_MAX_ITER 150000
#define FULL_IDLE_MS    900      /* idle time before the full render      */
#define SPIN_STEP       3.0      /* degrees of yaw per spin frame         */
#define UD_OFS          100      /* up-down id = its edit's id + UD_OFS   */

/* Program version, kept in sync with the chaoscope.c renderer. */
#define CHAOSCOPE_VERSION "1.0.0"

enum {
    IDC_SYS = 100, IDC_DESC,
    IDC_PLAB = 110,                       /* 110..115 parameter labels      */
    IDC_PED  = 120,                       /* 120..125 parameter edits       */
    IDC_DT   = 140, IDC_DEFAULTS,
    IDC_YAW  = 150, IDC_PITCH, IDC_ROLL, IDC_ZOOM, IDC_CYCLE, IDC_GAMMA,   /* sliders  */
    IDC_YAWV = 160,                       /* 160..164 value labels          */
    IDC_RESETVIEW = 170,
    IDC_PAL  = 180, IDC_SPLAT, IDC_W, IDC_H, IDC_ITER, IDC_SEED, IDC_DICE,
    IDC_SPIN, IDC_AUTO,
    IDC_RENDER = 190, IDC_SAVE, IDC_COPY, IDC_FRAMES, IDC_ANIM,
    IDC_STATUS = 300
};

enum { T_POLL = 1, T_FULL = 2 };

#define S(v)   MulDiv((v), g.dpi, 96)
#define CTL(id) GetDlgItem(g.panel, (id))

/* ------------------------------------------------------------------ */
/* Data                                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    char   name[32], desc[128];
    double dt;
    int    np;                  /* parameters shown in the UI (<= UIP)    */
    char   pn[UIP][16];
    double pd[UIP];
} Sys;

typedef struct {
    int    sys, pal, splat, w, h, iters, seed;
    double p[UIP], dt;
    double yaw, pitch, roll, zoom, cycle, panx, pany, gam;
} Settings;

static Sys  g_sys[MAXSYS];  static int g_nsys;
static char g_pal[MAXPAL][24]; static int g_npal;

static const struct { int id; const char *label; int lo, hi; double scale; const char *fmt; } VIEW[6] = {
    { IDC_YAW,   "Yaw",   -180, 180, 1.0,  "%.0f\xB0" },
    { IDC_PITCH, "Pitch", -180, 180, 1.0,  "%.0f\xB0" },
    { IDC_ROLL,  "Roll",  -180, 180, 1.0,  "%.0f\xB0" },
    { IDC_ZOOM,  "Zoom",    10, 2000, 0.01, "%.2fx"    },
    { IDC_CYCLE, "Cycle",    0, 100, 0.01, "%.2f"     },
    { IDC_GAMMA, "Gamma",   30, 400, 0.01, "%.2f"     },
};
static const double VIEW_DEFAULT[6] = { 180.0, 180.0, 0.0, 1.0, 0.0, 1.0 };

typedef enum { K_PREVIEW = 1, K_FULL, K_SAVE, K_ANIM } Kind;

static struct {
    HINSTANCE inst;
    HWND   hwnd, panel, view, status;
    HFONT  font;
    HBRUSH bg;
    HBITMAP bmp;
    int    dpi, content_h, panel_w;
    int    ui_lock, pending, spin, autofull, full_valid;
    int    drag, dx0, dy0;
    double    panx, pany;
    unsigned seq;
    double vw[6];                       /* yaw, pitch, roll, zoom, cycle   */
    char   exe[MAX_PATH], tmpdir[MAX_PATH], fullfile[MAX_PATH];
} g;

static struct {
    int    active;
    Kind   kind;
    HANDLE proc, pipe;
    DWORD  t0;
    char   out[MAX_PATH];
    char   line[600]; int linelen;
    char   err[400];
    int    frames_done, frames_total;
    Settings st;
} job;

static void request_render(void);
static void start_job(Kind k, const char *path, int frames);

/* ------------------------------------------------------------------ */
/* Small helpers                                                       */
/* ------------------------------------------------------------------ */

static void set_status(const char *fmt, ...) {
    char b[600];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(b, sizeof b, fmt, ap);
    va_end(ap);
    b[sizeof b - 1] = 0;
    SendMessage(g.status, SB_SETTEXT, 0, (LPARAM)b);
}

/* Split the status bar so its right part always shows the version. */
static void status_setup(void) {
    int parts[2];
    RECT rc;
    if (!g.status) return;
    GetClientRect(g.hwnd, &rc);
    parts[0] = rc.right - S(250);
    if (parts[0] < 0) parts[0] = 0;
    parts[1] = -1;                             /* extend to the right edge */
    SendMessage(g.status, SB_SETPARTS, 2, (LPARAM)parts);
    SendMessage(g.status, SB_SETTEXT, 1, (LPARAM)"Chaoscope by Shpati Koleka, v" CHAOSCOPE_VERSION);
}

static int finite_d(double v) { return v == v && v > -1e300 && v < 1e300; }

static double wrap180(double a) {
    a = fmod(a + 180.0, 360.0);
    if (a < 0.0) a += 360.0;
    return a - 180.0;
}

static double nice_step(double v) {          /* 1/10 of the leading digit */
    double a = fabs(v);
    if (a < 1e-12) return 0.1;
    return pow(10.0, floor(log10(a)) - 1.0);
}

static int get_text(int id, char *b, int cap) {
    GetWindowText(CTL(id), b, cap);
    return (int)strlen(b);
}

static int get_dbl(int id, double *v) {
    char b[64], *e, *c;
    get_text(id, b, sizeof b);
    for (c = b; *c; c++) if (*c == ',') *c = '.';
    *v = strtod(b, &e);
    if (e == b) return 0;
    while (*e == ' ') e++;
    return *e == 0 && finite_d(*v);
}

static int get_int(int id, int *v) {
    char b[64], *e;
    long x;
    get_text(id, b, sizeof b);
    x = strtol(b, &e, 10);
    if (e == b) return 0;
    while (*e == ' ') e++;
    if (*e != 0 || x < -2000000000L || x > 2000000000L) return 0;
    *v = (int)x;
    return 1;
}

static void set_edit(int id, double v) {      /* silent (no re-render)    */
    char b[48];
    snprintf(b, sizeof b, "%.10g", v);
    g.ui_lock++;
    SetWindowText(CTL(id), b);
    g.ui_lock--;
}

static void set_edit_notify(int id, double v) {   /* triggers EN_CHANGE   */
    char b[48];
    snprintf(b, sizeof b, "%.10g", v);
    SetWindowText(CTL(id), b);
}

/* ------------------------------------------------------------------ */
/* Child processes                                                     */
/* ------------------------------------------------------------------ */

static int spawn(const char *cmdline, HANDLE *proc, HANDLE *pipe_rd) {
    SECURITY_ATTRIBUTES sa;
    HANDLE rd = NULL, wr = NULL, nul;
    STARTUPINFO si;
    PROCESS_INFORMATION pi;
    char *cmd;
    BOOL ok;

    sa.nLength = sizeof sa; sa.lpSecurityDescriptor = NULL; sa.bInheritHandle = TRUE;
    if (!CreatePipe(&rd, &wr, &sa, 0)) return 0;
    SetHandleInformation(rd, HANDLE_FLAG_INHERIT, 0);
    nul = CreateFile("NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &sa,
                     OPEN_EXISTING, 0, NULL);
    memset(&si, 0, sizeof si);
    si.cb = sizeof si;
    si.dwFlags = STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_HIDE;
    si.hStdInput = nul; si.hStdOutput = wr; si.hStdError = wr;

    cmd = (char *)malloc(strlen(cmdline) + 1);
    if (!cmd) { CloseHandle(rd); CloseHandle(wr); if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul); return 0; }
    strcpy(cmd, cmdline);
    ok = CreateProcess(NULL, cmd, NULL, NULL, TRUE, CREATE_NO_WINDOW, NULL, NULL, &si, &pi);
    free(cmd);
    CloseHandle(wr);                       /* so the pipe reports EOF     */
    if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul);
    if (!ok) { CloseHandle(rd); return 0; }
    CloseHandle(pi.hThread);
    *proc = pi.hProcess; *pipe_rd = rd;
    return 1;
}

static DWORD pipe_read(HANDLE pipe, char *buf, DWORD cap) {
    DWORD avail = 0, got = 0;
    if (!PeekNamedPipe(pipe, NULL, 0, NULL, &avail, NULL) || avail == 0) return 0;
    if (avail > cap) avail = cap;
    if (!ReadFile(pipe, buf, avail, &got, NULL)) return 0;
    return got;
}

/* Run to completion, collecting stdout+stderr (used for --list-tsv). */
static int run_capture(const char *cmd, char *out, size_t cap, DWORD timeout_ms) {
    HANDLE proc, pipe;
    size_t n = 0;
    DWORD t0 = GetTickCount(), got;
    if (!spawn(cmd, &proc, &pipe)) return 0;
    for (;;) {
        got = pipe_read(pipe, out + n, (DWORD)(cap - 1 - n));
        n += got;
        if (n >= cap - 1) break;
        if (!got && WaitForSingleObject(proc, 10) == WAIT_OBJECT_0) {
            while ((got = pipe_read(pipe, out + n, (DWORD)(cap - 1 - n))) > 0) n += got;
            break;
        }
        if (GetTickCount() - t0 > timeout_ms) { TerminateProcess(proc, 1); break; }
    }
    out[n] = 0;
    CloseHandle(proc); CloseHandle(pipe);
    return n > 0;
}

static int split_tabs(char *line, char **f, int maxf) {
    int n = 0;
    char *c;
    f[n++] = line;
    for (c = line; *c && n < maxf; c++)
        if (*c == '\t') { *c = 0; f[n++] = c + 1; }
    return n;
}

static int load_systems(void) {
    static char buf[1 << 16];
    char cmd[MAX_PATH + 32], *line, *next;
    snprintf(cmd, sizeof cmd, "\"%s\" --list-tsv", g.exe);
    if (!run_capture(cmd, buf, sizeof buf, 8000)) return 0;
    for (line = buf; line && *line; line = next) {
        char *f[40];
        int nf, i;
        size_t L;
        next = strchr(line, '\n');
        if (next) *next++ = 0;
        L = strlen(line);
        if (L && line[L - 1] == '\r') line[L - 1] = 0;
        nf = split_tabs(line, f, 40);
        if (f[0][0] == 'S' && nf >= 5 && g_nsys < MAXSYS) {
            Sys *s = &g_sys[g_nsys];
            int np = atoi(f[4]);
            memset(s, 0, sizeof *s);
            lstrcpyn(s->name, f[1], sizeof s->name);
            lstrcpyn(s->desc, f[2], sizeof s->desc);
            s->dt = strtod(f[3], NULL);
            if (np > UIP) np = UIP;
            if (np > (nf - 5) / 2) np = (nf - 5) / 2;
            s->np = np < 0 ? 0 : np;
            for (i = 0; i < s->np; i++) {
                lstrcpyn(s->pn[i], f[5 + 2 * i], sizeof s->pn[i]);
                s->pd[i] = strtod(f[6 + 2 * i], NULL);
            }
            g_nsys++;
        } else if (f[0][0] == 'P' && nf >= 2 && g_npal < MAXPAL) {
            lstrcpyn(g_pal[g_npal++], f[1], sizeof g_pal[0]);
        }
    }
    return g_nsys > 0 && g_npal > 0;
}

/* ------------------------------------------------------------------ */
/* Settings <-> command line                                           */
/* ------------------------------------------------------------------ */

static int read_settings(Settings *s, char *err, size_t ecap) {
    const Sys *sy;
    int i, v;
    memset(s, 0, sizeof *s);
    s->sys = (int)SendMessage(CTL(IDC_SYS), CB_GETCURSEL, 0, 0);
    if (s->sys < 0 || s->sys >= g_nsys) { snprintf(err, ecap, "No system selected"); return 0; }
    sy = &g_sys[s->sys];
    for (i = 0; i < sy->np; i++)
        if (!get_dbl(IDC_PED + i, &s->p[i])) { snprintf(err, ecap, "Invalid number for parameter %s", sy->pn[i]); return 0; }
    if (!get_dbl(IDC_DT, &s->dt) || s->dt <= 0.0) { snprintf(err, ecap, "dt must be a number > 0"); return 0; }
    if (!get_int(IDC_W, &s->w) || s->w < 16 || s->w > 16384) { snprintf(err, ecap, "Width must be 16..16384"); return 0; }
    if (!get_int(IDC_H, &s->h) || s->h < 16 || s->h > 16384) { snprintf(err, ecap, "Height must be 16..16384"); return 0; }
    if (!get_int(IDC_ITER, &s->iters) || s->iters < 1000) { snprintf(err, ecap, "Iterations must be >= 1000"); return 0; }
    if (!get_int(IDC_SEED, &v) || v < 1) { snprintf(err, ecap, "Seed must be a positive integer"); return 0; }
    s->seed = v;
    s->splat = (int)SendMessage(CTL(IDC_SPLAT), CB_GETCURSEL, 0, 0);
    if (s->splat < 0) s->splat = 2;
    s->pal = (int)SendMessage(CTL(IDC_PAL), CB_GETCURSEL, 0, 0);
    if (s->pal < 0) s->pal = 0;
    s->yaw = g.vw[0]; s->pitch = g.vw[1]; s->roll = g.vw[2];
    s->zoom = g.vw[3]; s->cycle = g.vw[4]; s->gam = g.vw[5];
    s->panx = g.panx; s->pany = g.pany;
    return 1;
}

static void to_preview(Settings *s) {
    int m = s->w > s->h ? s->w : s->h;
    double sc = m > PREVIEW_MAX_DIM ? (double)PREVIEW_MAX_DIM / m : 1.0;
    s->w = (int)(s->w * sc + 0.5); if (s->w < 16) s->w = 16;
    s->h = (int)(s->h * sc + 0.5); if (s->h < 16) s->h = 16;
    s->splat = (int)(s->splat * sc + 0.5);
    s->panx *= sc;
    s->pany *= sc;
    if (s->iters > PREVIEW_MAX_ITER) s->iters = PREVIEW_MAX_ITER;
}

static void build_cmd(char *d, size_t cap, const Settings *s, const char *exe,
                      const char *out, int frames, const char *prefix) {
    const Sys *sy = &g_sys[s->sys];
    int n, i;
    n = snprintf(d, cap,
        "\"%s\" --system %s --width %d --height %d --iterations %d --dt %.10g"
        " --yaw %.6g --pitch %.6g --roll %.6g --zoom %.6g --palette %s"
        " --cycle %.6g --splat %d --seed %d --panx %.6g --pany %.6g --gamma %.4g",
        exe, sy->name, s->w, s->h, s->iters, s->dt, s->yaw, s->pitch, s->roll,
        s->zoom, g_pal[s->pal], s->cycle, s->splat, s->seed, s->panx, s->pany, s->gam);
    if (n < 0 || (size_t)n >= cap) return;
    if (sy->np > 0) {
        n += snprintf(d + n, cap - n, " --params ");
        for (i = 0; i < sy->np && (size_t)n < cap; i++)
            n += snprintf(d + n, cap - n, "%s%.10g", i ? "," : "", s->p[i]);
    }
    if ((size_t)n >= cap) return;
    if (frames > 0)
        snprintf(d + n, cap - n, " --animation %d --prefix \"%s\"", frames, prefix);
    else
        snprintf(d + n, cap - n, " --output \"%s\"", out);
}

/* ------------------------------------------------------------------ */
/* Render jobs                                                         */
/* ------------------------------------------------------------------ */

static int show_file(const char *path) {
    HBITMAP b = (HBITMAP)LoadImage(NULL, path, IMAGE_BITMAP, 0, 0,
                                   LR_LOADFROMFILE | LR_CREATEDIBSECTION);
    if (!b) return 0;
    if (g.bmp) DeleteObject(g.bmp);
    g.bmp = b;
    InvalidateRect(g.view, NULL, FALSE);
    return 1;
}

static void anim_button_idle(void) {
    SetWindowText(CTL(IDC_ANIM), "Render animation...");
}

static void job_kill(void) {
    if (!job.active) return;
    TerminateProcess(job.proc, 1);
    WaitForSingleObject(job.proc, 500);
    CloseHandle(job.proc); CloseHandle(job.pipe);
    job.active = 0;
    KillTimer(g.hwnd, T_POLL);
    if (job.kind == K_PREVIEW || job.kind == K_FULL) DeleteFile(job.out);
    if (job.kind == K_ANIM) anim_button_idle();
}

static void job_line(const char *line) {
    if (strncmp(line, "wrote ", 6) == 0) {
        if (job.kind == K_ANIM) {
            job.frames_done++;
            show_file(line + 6);
            set_status("Animation: frame %d of %d", job.frames_done, job.frames_total);
        }
        return;
    }
    {
        size_t n = strlen(job.err), L = strlen(line);
        if (L && n + L + 2 < sizeof job.err) {
            if (n) job.err[n++] = ' ';
            memcpy(job.err + n, line, L + 1);
        }
    }
}

static void job_feed(const char *d, DWORD n) {
    DWORD i;
    for (i = 0; i < n; i++) {
        char ch = d[i];
        if (ch == '\n') { job.line[job.linelen] = 0; job_line(job.line); job.linelen = 0; }
        else if (ch != '\r' && job.linelen < (int)sizeof job.line - 1) job.line[job.linelen++] = ch;
    }
}

static void job_finish(void) {
    DWORD code = 1, got, el = GetTickCount() - job.t0;
    char buf[512], out[MAX_PATH];
    Kind k = job.kind;
    Settings st = job.st;

    GetExitCodeProcess(job.proc, &code);
    while ((got = pipe_read(job.pipe, buf, sizeof buf)) > 0) job_feed(buf, got);
    if (job.linelen) { job.line[job.linelen] = 0; job_line(job.line); job.linelen = 0; }
    CloseHandle(job.proc); CloseHandle(job.pipe);
    job.active = 0;
    KillTimer(g.hwnd, T_POLL);
    lstrcpyn(out, job.out, MAX_PATH);
    if (k == K_ANIM) anim_button_idle();

    if (code != 0) {
        set_status("%s", job.err[0] ? job.err : "chaoscope.exe failed");
        if (k == K_PREVIEW || k == K_FULL) DeleteFile(out);
        if (g.spin) { g.spin = 0; SendMessage(CTL(IDC_SPIN), BM_SETCHECK, BST_UNCHECKED, 0); }
    } else {
        const char *warn = job.err[0] ? job.err : NULL;
        switch (k) {
        case K_PREVIEW:
            if (!show_file(out)) set_status("Could not load the rendered image");
            else set_status("Preview %dx%d in %lu ms%s%s", st.w, st.h, (unsigned long)el,
                            warn ? "  -  " : "", warn ? warn : "");
            DeleteFile(out);
            break;
        case K_FULL:
            if (!show_file(out)) { set_status("Could not load the rendered image"); DeleteFile(out); break; }
            if (g.fullfile[0]) DeleteFile(g.fullfile);
            lstrcpyn(g.fullfile, out, MAX_PATH);
            g.full_valid = 1;
            set_status("%s: %dx%d, %d iterations, %.2f s%s%s", g_sys[st.sys].name, st.w, st.h,
                       st.iters, el / 1000.0, warn ? "  -  " : "", warn ? warn : "");
            break;
        case K_SAVE:
            show_file(out);
            set_status("Saved %s", out);
            break;
        case K_ANIM:
            set_status("Animation finished: %d frames in %.1f s", job.frames_done, el / 1000.0);
            break;
        }
    }

    if (k == K_PREVIEW) {                       /* chain follow-up previews */
        if (g.pending) { g.pending = 0; start_job(K_PREVIEW, NULL, 0); }
        else if (g.spin && code == 0) {
            g.vw[0] = wrap180(g.vw[0] + SPIN_STEP);
            {   int i = 0;
                SendMessage(CTL(VIEW[i].id), TBM_SETPOS, TRUE, (LPARAM)(int)floor(g.vw[i] / VIEW[i].scale + 0.5));
            }
            { char b[32]; snprintf(b, sizeof b, VIEW[0].fmt, g.vw[0]); SetWindowText(CTL(IDC_YAWV), b); }
            start_job(K_PREVIEW, NULL, 0);
        }
    }
}

static void job_poll(void) {
    char buf[512];
    DWORD got;
    if (!job.active) { KillTimer(g.hwnd, T_POLL); return; }
    while ((got = pipe_read(job.pipe, buf, sizeof buf)) > 0) job_feed(buf, got);
    if (WaitForSingleObject(job.proc, 0) == WAIT_OBJECT_0) job_finish();
}

static void start_job(Kind k, const char *path, int frames) {
    Settings s;
    char err[160], cmd[2048], out[MAX_PATH];
    if (!read_settings(&s, err, sizeof err)) { set_status("%s", err); return; }
    job_kill();
    g.seq++;
    out[0] = 0;
    if (k == K_PREVIEW) to_preview(&s);
    if (k == K_PREVIEW || k == K_FULL)
        snprintf(out, sizeof out, "%sccgui_%lu_%u.bmp", g.tmpdir, (unsigned long)GetCurrentProcessId(), g.seq);
    else if (k == K_SAVE && path)
        lstrcpyn(out, path, MAX_PATH);
    build_cmd(cmd, sizeof cmd, &s, g.exe, out, k == K_ANIM ? frames : 0, path ? path : "");

    memset(&job, 0, sizeof job);
    if (!spawn(cmd, &job.proc, &job.pipe)) { set_status("Cannot start %s", g.exe); return; }
    job.active = 1; job.kind = k; job.t0 = GetTickCount(); job.st = s;
    job.frames_total = frames;
    lstrcpyn(job.out, out, MAX_PATH);
    SetTimer(g.hwnd, T_POLL, 15, NULL);
    if (k == K_FULL) set_status("Rendering %s %dx%d, %d iterations...", g_sys[s.sys].name, s.w, s.h, s.iters);
    else if (k == K_SAVE) set_status("Rendering to %s ...", out);
    else if (k == K_ANIM) set_status("Animation: starting (%d frames)...", frames);
}

/* A control changed: quick preview now, full render when things settle. */
static void request_render(void) {
    if (g.ui_lock) return;
    if (job.active && job.kind == K_ANIM) return;
    g.full_valid = 0;
    if (g.autofull && !g.spin) SetTimer(g.hwnd, T_FULL, FULL_IDLE_MS, NULL);
    else KillTimer(g.hwnd, T_FULL);
    if (job.active && job.kind == K_PREVIEW) { g.pending = 1; return; }
    start_job(K_PREVIEW, NULL, 0);
}

static void render_full_now(void) {
    KillTimer(g.hwnd, T_FULL);
    if (job.active && job.kind == K_ANIM) { set_status("An animation is rendering - cancel it first"); return; }
    g.pending = 0;
    start_job(K_FULL, NULL, 0);
}

/* ------------------------------------------------------------------ */
/* UI state                                                            */
/* ------------------------------------------------------------------ */

static void view_label(int i) {
    char b[32];
    snprintf(b, sizeof b, VIEW[i].fmt, g.vw[i]);
    SetWindowText(CTL(IDC_YAWV + i), b);
}

static void set_view_ui(void) {
    int i;
    g.ui_lock++;
    for (i = 0; i < 6; i++) {
        SendMessage(CTL(VIEW[i].id), TBM_SETPOS, TRUE, (LPARAM)(int)floor(g.vw[i] / VIEW[i].scale + 0.5));
        view_label(i);
    }
    g.ui_lock--;
}

static void apply_system(int si) {
    const Sys *sy = &g_sys[si];
    char b[200];
    int i;
    g.ui_lock++;
    snprintf(b, sizeof b, "%s%s", sy->desc, sy->np == 0 ? "  (no free parameters)" : "");
    SetWindowText(CTL(IDC_DESC), b);
    for (i = 0; i < UIP; i++) {
        int show = i < sy->np;
        ShowWindow(CTL(IDC_PLAB + i), show ? SW_SHOW : SW_HIDE);
        ShowWindow(CTL(IDC_PED + i), show ? SW_SHOW : SW_HIDE);
        ShowWindow(CTL(IDC_PED + i + UD_OFS), show ? SW_SHOW : SW_HIDE);
        if (show) { SetWindowText(CTL(IDC_PLAB + i), sy->pn[i]); set_edit(IDC_PED + i, sy->pd[i]); }
    }
    set_edit(IDC_DT, sy->dt);
    g.ui_lock--;
}

/* ------------------------------------------------------------------ */
/* Scrolling control panel                                             */
/* ------------------------------------------------------------------ */

static void panel_update_scroll(void) {
    RECT rc;
    SCROLLINFO si;
    GetClientRect(g.panel, &rc);
    si.cbSize = sizeof si; si.fMask = SIF_RANGE | SIF_PAGE;
    si.nMin = 0; si.nMax = g.content_h - 1; si.nPage = (UINT)rc.bottom;
    SetScrollInfo(g.panel, SB_VERT, &si, TRUE);
}

static void panel_scroll_to(int pos) {
    SCROLLINFO si;
    int old, maxpos;
    si.cbSize = sizeof si; si.fMask = SIF_ALL;
    GetScrollInfo(g.panel, SB_VERT, &si);
    maxpos = si.nMax - (int)si.nPage + 1;
    if (maxpos < 0) maxpos = 0;
    if (pos > maxpos) pos = maxpos;
    if (pos < 0) pos = 0;
    old = si.nPos;
    if (pos == old) return;
    si.fMask = SIF_POS; si.nPos = pos;
    SetScrollInfo(g.panel, SB_VERT, &si, TRUE);
    ScrollWindowEx(g.panel, 0, old - pos, NULL, NULL, NULL, NULL,
                   SW_SCROLLCHILDREN | SW_INVALIDATE | SW_ERASE);
    UpdateWindow(g.panel);
}

static LRESULT CALLBACK PanelProc(HWND h, UINT m, WPARAM w, LPARAM l) {
    switch (m) {
    case WM_COMMAND: case WM_HSCROLL: case WM_NOTIFY:
        return SendMessage(GetParent(h), m, w, l);
    case WM_VSCROLL: {
        SCROLLINFO si;
        int pos;
        si.cbSize = sizeof si; si.fMask = SIF_ALL;
        GetScrollInfo(h, SB_VERT, &si);
        pos = si.nPos;
        switch (LOWORD(w)) {
        case SB_LINEUP:    pos -= S(24); break;
        case SB_LINEDOWN:  pos += S(24); break;
        case SB_PAGEUP:    pos -= (int)si.nPage; break;
        case SB_PAGEDOWN:  pos += (int)si.nPage; break;
        case SB_TOP:       pos = 0; break;
        case SB_BOTTOM:    pos = si.nMax; break;
        case SB_THUMBTRACK: case SB_THUMBPOSITION: pos = si.nTrackPos; break;
        }
        panel_scroll_to(pos);
        return 0;
    }
    case WM_MOUSEWHEEL: {
        SCROLLINFO si;
        si.cbSize = sizeof si; si.fMask = SIF_POS;
        GetScrollInfo(h, SB_VERT, &si);
        panel_scroll_to(si.nPos - (short)HIWORD(w) * S(48) / 120);
        return 0;
    }
    }
    return DefWindowProc(h, m, w, l);
}

/* ------------------------------------------------------------------ */
/* Image view                                                          */
/* ------------------------------------------------------------------ */

static void paint_view(HWND h) {
    PAINTSTRUCT ps;
    HDC dc = BeginPaint(h, &ps);
    RECT rc;
    GetClientRect(h, &rc);
    if (g.bmp && rc.right > 0 && rc.bottom > 0) {
        BITMAP bm;
        GetObject(g.bmp, sizeof bm, &bm);
        if (bm.bmWidth > 0 && bm.bmHeight > 0) {
            double sx = (double)rc.right / bm.bmWidth, sy = (double)rc.bottom / bm.bmHeight;
            double sc = sx < sy ? sx : sy;
            int w = (int)(bm.bmWidth * sc), hh = (int)(bm.bmHeight * sc);
            int x, y;
            HDC mem;
            HGDIOBJ old;
            if (w < 1) w = 1;
            if (hh < 1) hh = 1;
            x = (rc.right - w) / 2; y = (rc.bottom - hh) / 2;
            mem = CreateCompatibleDC(dc);
            old = SelectObject(mem, g.bmp);
            SetStretchBltMode(dc, HALFTONE);
            SetBrushOrgEx(dc, 0, 0, NULL);
            StretchBlt(dc, x, y, w, hh, mem, 0, 0, bm.bmWidth, bm.bmHeight, SRCCOPY);
            SelectObject(mem, old);
            DeleteDC(mem);
            ExcludeClipRect(dc, x, y, x + w, y + hh);
        }
    }
    FillRect(dc, &rc, g.bg);
    if (!g.bmp) {
        SetBkMode(dc, TRANSPARENT);
        SetTextColor(dc, RGB(150, 150, 160));
        SelectObject(dc, g.font);
        DrawText(dc, "Rendering...", -1, &rc, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    }
    EndPaint(h, &ps);
}

static LRESULT CALLBACK ViewProc(HWND h, UINT m, WPARAM w, LPARAM l) {
    switch (m) {
    case WM_ERASEBKGND: return 1;
    case WM_PAINT: paint_view(h); return 0;
    case WM_SETCURSOR:
        SetCursor(LoadCursor(NULL, IDC_SIZEALL));
        return TRUE;
    case WM_LBUTTONDOWN: case WM_RBUTTONDOWN: case WM_MBUTTONDOWN:
        SetFocus(h);
        SetCapture(h);
        g.drag = (m == WM_LBUTTONDOWN) ? 1 : (m == WM_RBUTTONDOWN) ? 2 : 3;
        g.dx0 = (int)(short)LOWORD(l); g.dy0 = (int)(short)HIWORD(l);
        return 0;
    case WM_MOUSEMOVE:
        if (g.drag) {
            int x = (int)(short)LOWORD(l), y = (int)(short)HIWORD(l);
            double dx = x - g.dx0, dy = y - g.dy0;
            g.dx0 = x; g.dy0 = y;
            if (g.drag == 1) { g.vw[2] = wrap180(g.vw[2] + dx * 0.5); g.vw[1] = wrap180(g.vw[1] - dy * 0.5); }
            else if (g.drag == 2) { g.vw[0] = wrap180(g.vw[0] + dx * 0.5); }
            else { g.panx += dx; g.pany += dy; }
            set_view_ui();
            request_render();
        }
        return 0;
    case WM_LBUTTONUP: case WM_RBUTTONUP: case WM_MBUTTONUP:
        if (g.drag) { g.drag = 0; ReleaseCapture(); }
        return 0;
    case WM_MOUSEWHEEL: {
        double z = g.vw[3] * pow(1.1, (short)HIWORD(w) / 120.0);
        if (z < 0.1) z = 0.1;
        if (z > 20.0) z = 20.0;
        g.vw[3] = z;
        set_view_ui();
        request_render();
        return 0;
    }
    }
    return DefWindowProc(h, m, w, l);
}

/* ------------------------------------------------------------------ */
/* Building the controls                                               */
/* ------------------------------------------------------------------ */

static HWND mk(const char *cls, const char *text, DWORD style, DWORD ex,
               int x, int y, int w, int h, int id) {
    HWND c = CreateWindowEx(ex, cls, text, WS_CHILD | WS_VISIBLE | style,
                            S(x), S(y), S(w), S(h), g.panel, (HMENU)(INT_PTR)id, g.inst, NULL);
    SendMessage(c, WM_SETFONT, (WPARAM)g.font, TRUE);
    return c;
}

static void lab(const char *t, int x, int y, int w) { mk("STATIC", t, SS_LEFT, 0, x, y + 3, w, 18, 0); }
static void group(const char *t, int y, int h) { mk("BUTTON", t, BS_GROUPBOX | WS_CLIPSIBLINGS, WS_EX_TRANSPARENT, 6, y, 328, h, 0); }

static HWND edit(int id, int x, int y, int w) {
    return mk("EDIT", "", WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE, x, y, w, 22, id);
}

static void spinner(int edit_id) {
    HWND ud = CreateWindowEx(0, UPDOWN_CLASS, "",
                             WS_CHILD | WS_VISIBLE | UDS_ALIGNRIGHT | UDS_ARROWKEYS | UDS_NOTHOUSANDS,
                             0, 0, 0, 0, g.panel, (HMENU)(INT_PTR)(edit_id + UD_OFS), g.inst, NULL);
    SendMessage(ud, UDM_SETBUDDY, (WPARAM)CTL(edit_id), 0);
    SendMessage(ud, UDM_SETRANGE32, (WPARAM)-1000000000, 1000000000);
    SendMessage(ud, UDM_SETPOS32, 0, 0);
}

static void build_controls(void) {
    int i, y;
    HWND c;

    /* system */
    group("System", 4, 82);
    c = mk("COMBOBOX", "", WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST, 0, 14, 22, 312, 300, IDC_SYS);
    for (i = 0; i < g_nsys; i++) SendMessage(c, CB_ADDSTRING, 0, (LPARAM)g_sys[i].name);
    SendMessage(c, CB_SETDROPPEDWIDTH, S(200), 0);
    mk("STATIC", "", SS_LEFT, 0, 14, 52, 312, 32, IDC_DESC);

    /* parameters */
    group("Parameters", 92, 136);
    for (i = 0; i < UIP; i++) {
        int x = 14 + (i % 2) * 158, yy = 112 + (i / 2) * 28;
        mk("STATIC", "", SS_LEFT, 0, x, yy + 3, 46, 18, IDC_PLAB + i);
        edit(IDC_PED + i, x + 46, yy, 100);
        spinner(IDC_PED + i);
    }
    lab("dt", 14, 196, 46);
    edit(IDC_DT, 60, 196, 100);
    spinner(IDC_DT);
    mk("BUTTON", "Defaults", BS_PUSHBUTTON | WS_TABSTOP, 0, 172, 196, 154, 24, IDC_DEFAULTS);

    /* view */
    group("View", 234, 224);
    for (i = 0; i < 6; i++) {
        int yy = 254 + i * 28;
        lab(VIEW[i].label, 14, yy, 44);
        c = mk(TRACKBAR_CLASS, "", WS_TABSTOP | TBS_HORZ | TBS_NOTICKS, 0, 58, yy, 200, 24, VIEW[i].id);
        SendMessage(c, TBM_SETRANGE, TRUE, MAKELONG(VIEW[i].lo, VIEW[i].hi));
        SendMessage(c, TBM_SETPAGESIZE, 0, i < 3 ? 15 : 10);
        mk("STATIC", "", SS_LEFT, 0, 264, yy + 3, 62, 18, IDC_YAWV + i);
    }
    mk("BUTTON", "Reset view", BS_PUSHBUTTON | WS_TABSTOP, 0, 14, 424, 100, 26, IDC_RESETVIEW);
    mk("STATIC", "drag: roll/pitch  right-drag: yaw", SS_LEFT, 0, 124, 429, 202, 18, 0);

    /* render */
    group("Render", 464, 136);
    lab("Palette", 14, 484, 46);
    c = mk("COMBOBOX", "", WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST, 0, 60, 484, 112, 200, IDC_PAL);
    for (i = 0; i < g_npal; i++) SendMessage(c, CB_ADDSTRING, 0, (LPARAM)g_pal[i]);
    lab("Splat", 184, 484, 34);
    c = mk("COMBOBOX", "", WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST, 0, 220, 484, 106, 200, IDC_SPLAT);
    for (i = 0; i <= 4; i++) { char b[4]; snprintf(b, sizeof b, "%d", i); SendMessage(c, CB_ADDSTRING, 0, (LPARAM)b); }
    lab("Size", 14, 512, 46);
    edit(IDC_W, 60, 512, 56);
    lab("x", 121, 512, 10);
    edit(IDC_H, 134, 512, 56);
    lab("Iter", 200, 512, 30);
    edit(IDC_ITER, 232, 512, 94);
    lab("Seed", 14, 540, 46);
    edit(IDC_SEED, 60, 540, 70);
    mk("BUTTON", "New", BS_PUSHBUTTON | WS_TABSTOP, 0, 136, 539, 50, 24, IDC_DICE);
    mk("BUTTON", "Spin", BS_AUTOCHECKBOX | WS_TABSTOP, 0, 210, 541, 116, 20, IDC_SPIN);
    mk("BUTTON", "Full quality when idle", BS_AUTOCHECKBOX | WS_TABSTOP, 0, 14, 568, 312, 20, IDC_AUTO);

    /* actions */
    y = 608;
    mk("BUTTON", "Render", BS_DEFPUSHBUTTON | WS_TABSTOP, 0, 6, y, 100, 30, IDC_RENDER);
    mk("BUTTON", "Save image...", BS_PUSHBUTTON | WS_TABSTOP, 0, 112, y, 112, 30, IDC_SAVE);
    mk("BUTTON", "Copy command", BS_PUSHBUTTON | WS_TABSTOP, 0, 230, y, 104, 30, IDC_COPY);
    y = 646;
    lab("Frames", 6, y + 1, 46);
    edit(IDC_FRAMES, 54, y + 1, 56);
    mk("BUTTON", "Render animation...", BS_PUSHBUTTON | WS_TABSTOP, 0, 118, y, 216, 28, IDC_ANIM);
    g.content_h = S(y + 28 + 10);

    /* initial values */
    for (i = 0; i < 6; i++) g.vw[i] = VIEW_DEFAULT[i];
    set_view_ui();
    SendMessage(CTL(IDC_SPLAT), CB_SETCURSEL, 2, 0);
    SetWindowText(CTL(IDC_W), "3840");
    SetWindowText(CTL(IDC_H), "2160");
    SetWindowText(CTL(IDC_ITER), "500000");
    SetWindowText(CTL(IDC_SEED), "42");
    SetWindowText(CTL(IDC_FRAMES), "60");
    g.autofull = 1;
    SendMessage(CTL(IDC_AUTO), BM_SETCHECK, BST_CHECKED, 0);
    {
        int pi = 0, si = 0;
        for (i = 0; i < g_npal; i++) if (strcmp(g_pal[i], "chaoscope") == 0) pi = i;
        for (i = 0; i < g_nsys; i++) if (strcmp(g_sys[i].name, "lorenz") == 0) si = i;
        SendMessage(CTL(IDC_PAL), CB_SETCURSEL, pi, 0);
        SendMessage(CTL(IDC_SYS), CB_SETCURSEL, si, 0);
        apply_system(si);
    }
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static void copy_text(const char *t) {
    size_t n = strlen(t) + 1;
    HGLOBAL hg = GlobalAlloc(GMEM_MOVEABLE, n);
    if (!hg) return;
    memcpy(GlobalLock(hg), t, n);
    GlobalUnlock(hg);
    if (OpenClipboard(g.hwnd)) {
        EmptyClipboard();
        if (!SetClipboardData(CF_TEXT, hg)) GlobalFree(hg);
        CloseClipboard();
    } else GlobalFree(hg);
}

static int pick_save_path(char *path, const char *title, const char *def) {
    OPENFILENAME ofn;
    lstrcpyn(path, def, MAX_PATH);
    memset(&ofn, 0, sizeof ofn);
    ofn.lStructSize = sizeof ofn;
    ofn.hwndOwner = g.hwnd;
    ofn.lpstrFilter = "Bitmap (*.bmp)\0*.bmp\0\0";
    ofn.lpstrFile = path; ofn.nMaxFile = MAX_PATH;
    ofn.lpstrTitle = title;
    ofn.lpstrDefExt = "bmp";
    ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY;
    return GetSaveFileName(&ofn);
}

static void on_save(void) {
    char path[MAX_PATH], def[64];
    int si = (int)SendMessage(CTL(IDC_SYS), CB_GETCURSEL, 0, 0);
    snprintf(def, sizeof def, "%s.bmp", si >= 0 ? g_sys[si].name : "attractor");
    if (!pick_save_path(path, "Save image", def)) return;
    if (g.full_valid && g.fullfile[0]) {
        if (CopyFile(g.fullfile, path, FALSE)) set_status("Saved %s", path);
        else set_status("Could not write %s", path);
    } else {
        KillTimer(g.hwnd, T_FULL);
        if (job.active && job.kind == K_ANIM) { set_status("An animation is rendering - cancel it first"); return; }
        start_job(K_SAVE, path, 0);
    }
}

static void on_anim(void) {
    char path[MAX_PATH], *dot;
    int n;
    if (job.active && job.kind == K_ANIM) { job_kill(); set_status("Animation cancelled after %d frames", job.frames_done); return; }
    if (!get_int(IDC_FRAMES, &n) || n < 2 || n > 3600) { set_status("Frames must be 2..3600"); return; }
    if (!pick_save_path(path, "Frame file name (frames are saved as <name>_0000.bmp, ...)", "frame.bmp")) return;
    dot = strrchr(path, '.');
    if (dot && dot > strrchr(path, '\\')) *dot = 0;     /* chaoscope adds _NNNN.bmp */
    KillTimer(g.hwnd, T_FULL);
    start_job(K_ANIM, path, n);
    if (job.active) SetWindowText(CTL(IDC_ANIM), "Cancel animation");
}

static void on_copy(void) {
    Settings s;
    char err[160], cmd[2048];
    if (!read_settings(&s, err, sizeof err)) { set_status("%s", err); return; }
    build_cmd(cmd, sizeof cmd, &s, "chaoscope.exe", "attractor.bmp", 0, NULL);
    copy_text(cmd);
    set_status("Command line copied to the clipboard");
}

static int is_setting_edit(int id) {
    return (id >= IDC_PED && id < IDC_PED + UIP) || id == IDC_DT || id == IDC_W ||
           id == IDC_H || id == IDC_ITER || id == IDC_SEED;
}

static void on_command(WPARAM w) {
    int id = LOWORD(w), code = HIWORD(w);
    if (id == IDOK) { render_full_now(); return; }
    if (id == IDCANCEL || g.ui_lock) return;
    switch (id) {
    case IDC_SYS:
        if (code == CBN_SELCHANGE) {
            apply_system((int)SendMessage(CTL(IDC_SYS), CB_GETCURSEL, 0, 0));
            request_render();
        }
        return;
    case IDC_PAL: case IDC_SPLAT:
        if (code == CBN_SELCHANGE) request_render();
        return;
    case IDC_DEFAULTS:
        apply_system((int)SendMessage(CTL(IDC_SYS), CB_GETCURSEL, 0, 0));
        request_render();
        return;
    case IDC_RESETVIEW: {
        int i;
        for (i = 0; i < 6; i++) g.vw[i] = VIEW_DEFAULT[i];
        g.panx = g.pany = 0.0;
        set_view_ui();
        request_render();
        return;
    }
    case IDC_DICE:
        set_edit_notify(IDC_SEED, 1 + rand() % 999999);
        return;
    case IDC_SPIN:
        g.spin = SendMessage(CTL(IDC_SPIN), BM_GETCHECK, 0, 0) == BST_CHECKED;
        if (g.spin) { KillTimer(g.hwnd, T_FULL); if (!(job.active && job.kind == K_ANIM) && !(job.active && job.kind == K_PREVIEW)) start_job(K_PREVIEW, NULL, 0); }
        else request_render();
        return;
    case IDC_AUTO:
        g.autofull = SendMessage(CTL(IDC_AUTO), BM_GETCHECK, 0, 0) == BST_CHECKED;
        if (g.autofull) request_render(); else KillTimer(g.hwnd, T_FULL);
        return;
    case IDC_RENDER: render_full_now(); return;
    case IDC_SAVE:   on_save(); return;
    case IDC_COPY:   on_copy(); return;
    case IDC_ANIM:   on_anim(); return;
    }
    if (code == EN_CHANGE && is_setting_edit(id)) request_render();
}

static void on_hscroll(LPARAM l) {
    HWND tb = (HWND)l;
    int id = GetDlgCtrlID(tb), i;
    if (g.ui_lock) return;
    for (i = 0; i < 6; i++) {
        if (VIEW[i].id == id) {
            g.vw[i] = (int)SendMessage(tb, TBM_GETPOS, 0, 0) * VIEW[i].scale;
            view_label(i);
            request_render();
            return;
        }
    }
}

static LRESULT on_notify(NMHDR *nm) {
    if (nm->code == UDN_DELTAPOS) {
        NMUPDOWN *ud = (NMUPDOWN *)nm;
        int eid = (int)nm->idFrom - UD_OFS;
        double v, step;
        if (!get_dbl(eid, &v)) return 1;               /* unparsable: ignore */
        step = nice_step(v);
        if (GetKeyState(VK_SHIFT) & 0x8000) step *= 10.0;
        if (GetKeyState(VK_CONTROL) & 0x8000) step /= 10.0;
        v += (ud->iDelta > 0 ? 1.0 : -1.0) * step;
        if (eid == IDC_DT && v <= 0.0) return 1;
        set_edit_notify(eid, v);
        return 1;                                      /* we handled it     */
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Main window                                                         */
/* ------------------------------------------------------------------ */

static void layout(void) {
    RECT rc;
    int w, h, sh, pw;
    if (!g.panel || !g.status) return;
    GetClientRect(g.hwnd, &rc);
    w = rc.right; h = rc.bottom;
    SendMessage(g.status, WM_SIZE, 0, 0);
    status_setup();
    GetWindowRect(g.status, &rc);
    sh = rc.bottom - rc.top;
    pw = g.panel_w;
    if (pw > w) pw = w;
    MoveWindow(g.panel, 0, 0, pw, h - sh > 0 ? h - sh : 0, TRUE);
    MoveWindow(g.view, pw, 0, w - pw > 0 ? w - pw : 0, h - sh > 0 ? h - sh : 0, TRUE);
    InvalidateRect(g.view, NULL, FALSE); 
    panel_update_scroll();
}

static void cleanup_temp(void) {
    WIN32_FIND_DATA fd;
    char pat[MAX_PATH + 40], p[MAX_PATH + 40];
    HANDLE h;
    snprintf(pat, sizeof pat, "%sccgui_%lu_*.bmp", g.tmpdir, (unsigned long)GetCurrentProcessId());
    h = FindFirstFile(pat, &fd);
    if (h == INVALID_HANDLE_VALUE) return;
    do {
        snprintf(p, sizeof p, "%s%s", g.tmpdir, fd.cFileName);
        DeleteFile(p);
    } while (FindNextFile(h, &fd));
    FindClose(h);
}

static LRESULT CALLBACK MainProc(HWND h, UINT m, WPARAM w, LPARAM l) {
    switch (m) {
    case WM_CREATE:
        g.hwnd = h;
        g.status = CreateWindowEx(0, STATUSCLASSNAME, "", WS_CHILD | WS_VISIBLE | SBARS_SIZEGRIP,
                                  0, 0, 0, 0, h, (HMENU)(INT_PTR)IDC_STATUS, g.inst, NULL);
        SendMessage(g.status, WM_SETFONT, (WPARAM)g.font, TRUE);
        g.panel_w = S(340) + GetSystemMetrics(SM_CXVSCROLL);
        g.panel = CreateWindowEx(WS_EX_CONTROLPARENT, "CCGuiPanel", "",
                                 WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_CLIPCHILDREN,
                                 0, 0, g.panel_w, 100, h, NULL, g.inst, NULL);
        g.view = CreateWindowEx(0, "CCGuiView", "", WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS,
                                0, 0, 100, 100, h, NULL, g.inst, NULL);
        g.ui_lock++;
        build_controls();
        g.ui_lock--;
        layout();
        set_status("Drag the image to roll/pitch, right-drag to yaw, mouse wheel to zoom");
        return 0;
    case WM_SIZE:
        layout();
        return 0;
    case WM_GETMINMAXINFO:
        ((MINMAXINFO *)l)->ptMinTrackSize.x = S(760);
        ((MINMAXINFO *)l)->ptMinTrackSize.y = S(480);
        return 0;
    case WM_COMMAND:  on_command(w); return 0;
    case WM_HSCROLL:  on_hscroll(l); return 0;
    case WM_NOTIFY:   return on_notify((NMHDR *)l);
    case WM_TIMER:
        if (w == T_POLL) job_poll();
        else if (w == T_FULL) { KillTimer(h, T_FULL); if (!g.spin) render_full_now(); }
        return 0;
    case WM_DESTROY:
        job_kill();
        if (g.bmp) DeleteObject(g.bmp);
        cleanup_temp();
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProc(h, m, w, l);
}

static int locate_exe(const char *arg) {
    char *c;
    if (arg && arg[0]) {                         /* optional command-line path */
        lstrcpyn(g.exe, arg, MAX_PATH);
        if (g.exe[0] == '"') memmove(g.exe, g.exe + 1, strlen(g.exe));
        c = strrchr(g.exe, '"'); if (c) *c = 0;
        if (GetFileAttributes(g.exe) != INVALID_FILE_ATTRIBUTES) return 1;
    }
    GetModuleFileName(NULL, g.exe, MAX_PATH);
    c = strrchr(g.exe, '\\');
    if (c) c[1] = 0;
    lstrcat(g.exe, "chaoscope.exe");
    if (GetFileAttributes(g.exe) != INVALID_FILE_ATTRIBUTES) return 1;
    {
        char found[MAX_PATH];
        if (SearchPath(NULL, "chaoscope.exe", NULL, MAX_PATH, found, NULL)) { lstrcpyn(g.exe, found, MAX_PATH); return 1; }
    }
    {
        OPENFILENAME ofn;
        char path[MAX_PATH] = "";
        MessageBox(NULL, "chaoscope.exe was not found next to this program or on the PATH.\n"
                         "Please locate it.", "Chaoscope GUI", MB_OK | MB_ICONINFORMATION);
        memset(&ofn, 0, sizeof ofn);
        ofn.lStructSize = sizeof ofn;
        ofn.lpstrFilter = "chaoscope.exe\0chaoscope.exe\0All files\0*.*\0\0";
        ofn.lpstrFile = path; ofn.nMaxFile = MAX_PATH;
        ofn.lpstrTitle = "Locate chaoscope.exe";
        ofn.Flags = OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
        if (!GetOpenFileName(&ofn)) return 0;
        lstrcpyn(g.exe, path, MAX_PATH);
        return 1;
    }
}

int WINAPI WinMain(HINSTANCE inst, HINSTANCE prev, LPSTR cmd, int show) {
    INITCOMMONCONTROLSEX icc;
    WNDCLASS wc;
    HDC dc;
    HWND wnd;
    MSG msg;
    RECT work;
    int ww, wh;
    typedef BOOL (WINAPI *DPIFN)(void);
    DPIFN dpifn;
    (void)prev;

    dpifn = (DPIFN)(void *)GetProcAddress(GetModuleHandle("user32.dll"), "SetProcessDPIAware");
    if (dpifn) dpifn();                         /* no-op if the manifest did it */
    g.inst = inst;
    dc = GetDC(NULL); g.dpi = GetDeviceCaps(dc, LOGPIXELSX); ReleaseDC(NULL, dc);
    if (g.dpi < 96) g.dpi = 96;

    icc.dwSize = sizeof icc;
    icc.dwICC = ICC_BAR_CLASSES | ICC_UPDOWN_CLASS | ICC_STANDARD_CLASSES;
    InitCommonControlsEx(&icc);
    srand((unsigned)GetTickCount());

    GetTempPath(MAX_PATH, g.tmpdir);
    if (!locate_exe(cmd)) return 1;
    if (!load_systems()) {
        MessageBox(NULL, "Could not read the system list from chaoscope.exe --list-tsv.\n"
                         "Make sure chaoscope.exe was built from the matching chaoscope.c.",
                   "Chaoscope GUI", MB_OK | MB_ICONERROR);
        return 1;
    }

    g.font = CreateFont(-MulDiv(9, g.dpi, 72), 0, 0, 0, FW_NORMAL, 0, 0, 0, DEFAULT_CHARSET,
                        OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                        DEFAULT_PITCH | FF_DONTCARE, "Segoe UI");
    g.bg = CreateSolidBrush(RGB(0, 0, 0));

    memset(&wc, 0, sizeof wc);
    wc.hInstance = inst; wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_BTNFACE + 1);
    wc.lpszClassName = "CCGuiMain"; wc.lpfnWndProc = MainProc;
    wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    RegisterClass(&wc);
    wc.lpszClassName = "CCGuiPanel"; wc.lpfnWndProc = PanelProc; wc.hIcon = NULL;
    RegisterClass(&wc);
    wc.lpszClassName = "CCGuiView"; wc.lpfnWndProc = ViewProc; wc.hbrBackground = NULL;
    RegisterClass(&wc);

    SystemParametersInfo(SPI_GETWORKAREA, 0, &work, 0);
    ww = S(1200); wh = S(740);
    if (ww > work.right - work.left) ww = work.right - work.left;
    if (wh > work.bottom - work.top) wh = work.bottom - work.top;
    wnd = CreateWindowEx(0, "CCGuiMain", "Chaoscope", WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,
                         CW_USEDEFAULT, CW_USEDEFAULT, ww, wh, NULL, NULL, inst, NULL);
    if (!wnd) return 1;
    ShowWindow(wnd, show);
    UpdateWindow(wnd);
    request_render();

    while (GetMessage(&msg, NULL, 0, 0) > 0) {
        if (!IsDialogMessage(g.hwnd, &msg)) {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }
    return (int)msg.wParam;
}