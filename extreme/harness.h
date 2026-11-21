/* harness.h — 极端场景矩阵共享框架。
   每个用例在独立子进程中运行：崩溃（信号、sanitizer 中止、10s 挂起）
   被隔离并计数，不会中断整个矩阵。输出经 EXPECT_FINITE 检查 NaN/Inf。 */
#ifndef EXTREME_HARNESS_H
#define EXTREME_HARNESS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

static int g_bad; /* 当前用例内的违规计数（NaN / Inf / 断言失败） */

#define EXPECT_FINITE(name, expr) do {                                   \
    double v_ = (double)(expr);                                          \
    if (isnan(v_))      { printf("  NAN  %s\n", name); g_bad++; }        \
    else if (isinf(v_)) { printf("  INF  %s\n", name); g_bad++; }        \
} while (0)

#define EXPECT_TRUE(name, cond) do {                                     \
    if (!(cond)) { printf("  FAIL %s\n", name); g_bad++; }               \
} while (0)

typedef void (*case_fn)(void);
typedef struct { const char* category; const char* name; case_fn fn; } Case;

static int run_matrix(const char* module, const Case* cases, int ncases) {
    int crashes = 0, bad = 0;
    printf("module=%s cases=%d\n", module, ncases);
    for (int i = 0; i < ncases; i++) {
        g_bad = 0;
        fflush(stdout);
        pid_t pid = fork();
        if (pid == 0) {
            setvbuf(stdout, NULL, _IONBF, 0);
            alarm(10); /* 挂起 -> SIGALRM -> 记为崩溃 */
            cases[i].fn();
            _exit(g_bad == 0 ? 0 : (g_bad > 100 ? 100 : g_bad));
        }
        int status = 0, crashed = 0, nbad = 0;
        waitpid(pid, &status, 0);
        if (WIFSIGNALED(status)) {
            crashed = 1;
            printf("  SIGNAL %d\n", WTERMSIG(status));
        } else if (WIFEXITED(status)) {
            nbad = WEXITSTATUS(status);
        } else {
            crashed = 1;
        }
        crashes += crashed;
        bad += nbad;
        printf("%s/%s crash=%d violations=%d\n",
               cases[i].category, cases[i].name, crashed, nbad);
    }
    printf("summary module=%s crashes=%d violations=%d all_ok=%d\n",
           module, crashes, bad, crashes == 0 && bad == 0);
    return (crashes == 0 && bad == 0) ? 0 : 1;
}

#endif
