// Timing thresholds in tests are scaled by PERF_SLACK (set > 1 on slower CI runners).
const env = (globalThis as { process?: { env?: Record<string, string | undefined> } }).process?.env ?? {};
export const PERF_SLACK = Math.max(1, Number(env.PERF_SLACK ?? 1) || 1);
export const ms = (budget: number): number => budget * PERF_SLACK;
