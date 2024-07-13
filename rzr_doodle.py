#!/usr/bin/env python

import sys, os, time, argparse, subprocess, signal, base64, shutil

IMAGE_PATH = "_rzr_doodle.png"
EXECUTABLE = "_rzr_doodle"

class SignalRefreshViewer:
	def __init__(self, exe, sig):
		self.exe = exe
		self.path = shutil.which(exe)
		self.available = self.path is not None
		self.image_path = IMAGE_PATH
		self.process = None
		self.sig = sig

	def load(self):
		if self.process is None:
			if not self.available:
				sys.stderr.write("%s: not available\n" % self.exe)
				sys.exit(1)
			self.process = subprocess.Popen([self.path, self.image_path], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
		else:
			self.process.send_signal(self.sig)

	def did_exit(self):
		return self.process is not None and self.process.poll() is not None

viewers = []
viewers_not_found = []
def add_viewer(name, viewer):
	if not viewer.available:
		viewers_not_found.append(name)
		return
	viewers.append((name, viewer))
add_viewer("feh", SignalRefreshViewer("feh", signal.SIGUSR1))

if len(viewers) == 0:
	sys.stderr.write("No supported viewers are available (tried: %s)\n" % ",".join(viewers_not_found))
	sys.exit(1)


RED_DOT_PNG = base64.b64decode('''iVBORw0KGgoAAAANSUhEUgAAAQAAAAEAAQMAAABmvDolAAAAIGNIUk0AAHomAACAhAAA+gAAAIDoAAB1MAAA6mAAADqYAAAXcJy6UTwAAAAGUExURf8AAP///0EdNBEAAAABYktHRAH/Ai3eAAAAH0lEQVRo3u3BAQ0AAADCoPdPbQ43oAAAAAAAAAAAvg0hAAABfxmcpwAAAABJRU5ErkJggg==''') # red.png 256x256 red image

parser = argparse.ArgumentParser(prog = 'rzr_doodle')
parser.add_argument("c_source_files", nargs='+')
parser.add_argument("-W", "--width", type=int, help="default is 256")
parser.add_argument("-H", "--height", type=int, help="default is to use the width value")
parser.add_argument("-P", "--pixels-per-unit", type=float, help="default is half of width")
parser.add_argument("-S", "--supersample-factor", type=int, help="if N is passed, each pixel will be based on NxN samples; default is 12")
parser.add_argument("-I", "--ignore-unresolved-symbols", action='store_true', help="tells linker to ignore problems with unresolved symbols; useful if your source file has cascading dependencies, but also dangerous because /any/ access to an unresolved symbol crashes. you might prefer using `#ifndef RZR_DOODLE` blocks instead.")
parser.add_argument("-V", "--viewer", default="feh", type=str, help="available: %s" % ",".join([v[0] for v in viewers]))
args = parser.parse_args()

viewer = None
for v in viewers:
	if v[0] == args.viewer:
		viewer = v[1]
		break

if viewer is None:
	sys.stderr.write("viewer not available: %s\n", args.viewer)
	sys.exit(1)

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

DEFINES += [define("RZR_OUTPUT_IMAGE_PATH", '"'+IMAGE_PATH+'"')]


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
while not viewer.did_exit():
	builder.begin()
	builder.compile(ARTIFACT_MAIN, ["rzr.c", "rzr.h"], ["cc", "-O2"] + DEFINES + ["-c", "rzr.c", "-o", ARTIFACT_MAIN])
	for src in args.c_source_files:
		# "foo.c" => "_rzr_doodle_fn_foo.o"
		obj = "_rzr_doodle_fn_%s.o" % os.path.splitext(src)[0]
		builder.compile(obj, [src] + ["rzr.h"], ["cc", "-O0"] + DEFINES + ["-c", src, "-o", obj])
	builder.end()

