"""Abort the process if its working set exceeds a limit (Windows, ctypes; no psutil)."""
import ctypes
import os
import sys
import threading
import time
from ctypes import wintypes


class PMC(ctypes.Structure):
    _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]


def usage_mb():
    pmc = PMC(); pmc.cb = ctypes.sizeof(PMC)
    h = ctypes.windll.kernel32.GetCurrentProcess()
    ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.c_void_p(h), ctypes.byref(pmc), pmc.cb)
    return max(pmc.WorkingSetSize, pmc.PagefileUsage) / 2**20, pmc.PeakPagefileUsage / 2**20


state = {"label": ""}


def start(limit_mb=4600, on_abort=None):
    def run():
        while True:
            cur, _ = usage_mb()
            if cur > limit_mb:
                msg = f"MEMGUARD: {cur:.0f} MB > {limit_mb} MB during '{state['label']}', aborting"
                print(msg, flush=True)
                if on_abort:
                    on_abort(state["label"], cur)
                sys.stdout.flush()
                os._exit(3)
            time.sleep(0.3)
    threading.Thread(target=run, daemon=True).start()
