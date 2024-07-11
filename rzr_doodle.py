#!/usr/bin/env python

import sys, os, time, argparse, subprocess, signal, base64

RED_DOT_PNG = base64.b64decode('''iVBORw0KGgoAAAANSUhEUgAAAQAAAAEAAQMAAABmvDolAAAAIGNIUk0AAHomAACAhAAA+gAAAIDoAAB1MAAA6mAAADqYAAAXcJy6UTwAAAAGUExURf8AAP///0EdNBEAAAABYktHRAH/Ai3eAAAAH0lEQVRo3u3BAQ0AAADCoPdPbQ43oAAAAAAAAAAAvg0hAAABfxmcpwAAAABJRU5ErkJggg==''') # red.png 256x256 red image

parser = argparse.ArgumentParser(prog = 'rzr_doodle')
parser.add_argument("c_source_files", nargs='+')
parser.add_argument("-W", "--width", type=int)
parser.add_argument("-H", "--height", type=int)
parser.add_argument("-P", "--pixels-per-unit", type=float)
parser.add_argument("-S", "--supersample-factor", type=int)
parser.add_argument("-I", "--ignore-unresolved-symbols", action='store_true')
args = parser.parse_args()

def define(k, v=None):
	if v is None:
		return "-D%s" % k
	else:
		return "-D%s=%s" % (k,v)

DEFINES = [ define("RZR_DOODLE") ]
if args.width              is not None: DEFINES += [define("RZR_DOODLE_WIDTH", args.width)]
if args.height             is not None: DEFINES += [define("RZR_DOODLE_HEIGHT", args.height)]
if args.pixels_per_unit    is not None: DEFINES += [define("RZR_DOODLE_PIXELS_PER_UNIT", args.pixels_per_unit)]
if args.supersample_factor is not None: DEFINES += [define("RZR_SUPERSAMPLE_FACTOR", args.supersample_factor)]

IMAGE_PATH = "_rzr_doodle.png"
DEFINES += [define("RZR_OUTPUT_IMAGE_PATH", '"'+IMAGE_PATH+'"')]

class FEH:
	def __init__(self, path):
		self.path = path
		self.process = None

	def load(self):
		if self.process is None:
			self.process = subprocess.Popen(["feh", self.path])
		else:
			self.process.send_signal(signal.SIGUSR1)

viewer=FEH(IMAGE_PATH)

EXECUTABLE = "_rzr_doodle"

def mtime(path):
	try:
		return os.stat(path).st_mtime
	except FileNotFoundError as e:
		return 0

def run(argv):
	t0 = time.time()
	p = subprocess.run(argv)
	if p.returncode != 0:
		return False
	dt = time.time() - t0
	print("%s took %.3f seconds" % (argv, dt))
	return True

class Builder:
	def __init__(self):
		self.error_times = {}
		self.seen = set()
		self.showing_error = False

	def begin(self):
		self.objs = []
		self.error = False
		self.update = False

	def show_error(self):
		if self.showing_error: return
		with open(IMAGE_PATH,"wb") as f: f.write(RED_DOT_PNG)
		viewer.load()
		self.showing_error = True

	def wait(self, err=False):
		if err: self.show_error()
		time.sleep(0.05)

	def errwait(self):
		return self.wait(True)

	def end(self):
		if self.error: return self.errwait()
		if not self.update: return self.wait()
		extra = []
		if args.ignore_unresolved_symbols: extra += ["-Wl,--unresolved-symbols=ignore-in-object-files"]
		if not run(["cc"] + self.objs + ["-o", EXECUTABLE] + extra + ["-lm"]): return self.errwait()
		if not run(["./%s" % EXECUTABLE]): return self.errwait()
		if mtime(IMAGE_PATH) == 0: return self.errwait()
		viewer.load()
		self.showing_error = False

	def compile(self, artifact, dependencies, cc_argv):
		if self.error: return
		first = artifact not in self.seen
		self.objs += [artifact]
		ta = mtime(artifact)
		tdmax=0
		for d in dependencies:
			td = mtime(d)
			tdmax = max(tdmax,td)
		if not first and ta >= tdmax: return
		if artifact in self.error_times:
			et = self.error_times[artifact]
			if tdmax <= et:
				self.error = True
				return
			del self.error_times[artifact]
		if not run(cc_argv):
			self.error = True
			self.error_times[artifact] = tdmax
			return
		self.seen.add(artifact)
		self.update = True

builder = Builder()

ARTIFACT_MAIN = "_rzr_doodle_main.o"
while True:
	builder.begin()
	builder.compile(ARTIFACT_MAIN, ["rzr.c", "rzr.h"], ["cc", "-O2"] + DEFINES + ["-c", "rzr.c", "-o", ARTIFACT_MAIN])
	for src in args.c_source_files:
		# "foo.c" => "_rzr_doodle_fn_foo.o"
		obj = "_rzr_doodle_fn_%s.o" % os.path.splitext(src)[0]
		builder.compile(obj, [src] + ["rzr.h"], ["cc", "-O0"] + DEFINES + ["-c", src, "-o", obj])
	builder.end()

