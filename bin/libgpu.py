#!/usr/bin/env python

import os
import time
import subprocess as sp
import xml.etree.ElementTree as ET

class GPU:
    def __init__(self, host=None, login=None):
        if host is None:
            self.host = os.getenv('HOSTNAME')
            if self.host is None:
                self.host = sp.check_output("hostname").strip()
        else:
            self.host = host
        if login is None:
            self.login = self.host
        else:
            self.login = login

        self.n_gpu = 0
        self.proc_s = []
        self.check()
        self._visible_devices = None
    def __len__(self):
        return len(self.usable())
    def check(self):
        if self.host == None or self.host == self.login:
            self.check_proc = sp.Popen(["nvidia-smi", "-q", '-x'], stdout=sp.PIPE)
        else:
            self.check_proc = sp.Popen(["ssh", self.host, "nvidia-smi", "-q", "-x"], stdout=sp.PIPE)
    def update(self):
        t_wait0 = time.time()
        while True:
            if self.check_proc.poll() is not None:
                break
            if (time.time() - t_wait0) > 10.0:
                return
            time.sleep(0.1)
        xml = self.check_proc.stdout.read()
        info = ET.fromstring(xml)
        #
        self.n_gpu = int(info.find("attached_gpus").text)
        self.proc_s = []
        for gpu in info.findall("gpu"):
            proc_s = []
            for proc in gpu.find("processes"):
                proc_id = proc.find("pid").text
                proc_type = proc.find("type").text
                proc_name = proc.find("process_name").text.split()[0]
                if proc_type != 'G':
                    proc_s.append((proc_type, proc_id, proc_name))
            self.proc_s.append(proc_s)
        if self.host.startswith("yellow"):
            if len(self.proc_s) > 1:
                tmp = self.proc_s
                self.proc_s = [tmp[1], tmp[0]]
    @property
    def CUDA_VISIBLE_DEVICES(self):
        if self._visible_devices is not None:
            return self._visible_devices
        #
        if 'CUDA_VISIBLE_DEVICES' in os.environ:
            visible_devices = [int(x) for x in os.environ['CUDA_VISIBLE_DEVICES'].split(",")]
        else:
            visible_devices = list(range(len(self.proc_s)))
        self._visible_devices = visible_devices
        return visible_devices
    def usable(self, update=False):
        if update: self.update()
        #
        visible_devices = self.CUDA_VISIBLE_DEVICES
        #
        empty = []
        for k,proc_s in enumerate(self.proc_s):
            status = (k in visible_devices)
            for proc in proc_s:
                if proc[0] != 'G':  # graphic proc.
                    status = False
                    break
            if status:
                empty.append(k)
        if self.host.startswith("yellow"):
            empty.sort(reverse=True)
        return empty
