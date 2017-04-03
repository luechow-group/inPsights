#!/usr/bin/python

import sys

def read_configurations(file):
	norb = 0
	ne = 0
	nalpha = 0
	nbeta = 0
	configs = []
	with open(file, "r") as f:
		found = False
		for line in f:
			if "no. active" in line:
				lline = line.lower()
				if "orbitals" in lline:
					norb = int(line.split()[-1])
				elif "electrons" in lline:
					ne = int(line.split()[-1])
			elif "alpha electrons" in line:
				w = line.split()
				nalpha = int(w[0])
				nbeta = int(w[3])
			elif "Configuration" in line and "Symmetry" in line:
				found = True
				configs.append(line.split()[4])
			elif found:
				break

	if not norb or not ne or not configs:
		raise EOFError()

	assert nalpha + nbeta == ne

	virt_orbs = norb - (ne + 1) // 2
	assert virt_orbs >= 0

	dets = []
	for config in configs:
		offset = nalpha + virt_orbs - len(config) + 1
		det_alpha = range(1, offset)
		det_beta = det_alpha[:]

		for i, c in enumerate(config):
			if c == '1' or c == 'a':
				det_alpha.append(str(i + offset))
			if c == '1' or c == 'b':
				det_beta.append(str(i + offset))

		dets.append(' '.join(det_alpha + det_beta))

	return dets

def read_eigenvalues(file):
	eigenvalues = []
	with open(file, "r") as f:
		found = False
		for line in f:
			if line.startswith(" EIGENVECTOR USED TO COMPUTE DENSITY MATRICES"):
				eigenvalues = []
				found = True
				skip_next = True
			elif line.startswith(" Bra  VECTOR USED TO COMPUTE DENSITY MATRICES"):
				found = False
			elif found:
				skip_next = not skip_next
				if not skip_next: continue

				eigenvalues += map(float, line.replace("D", "E").split()[1:])

	return eigenvalues

def main(config_file, coeff_file):
	dets = read_configurations(config_file)
	print >> sys.stderr, "%d configurations found" % len(dets)
	eigenvalues = read_eigenvalues(coeff_file)
	print >> sys.stderr, "%d eigenvalues found" % len(eigenvalues)

	assert len(eigenvalues) >= len(dets)

	for ev, det in sorted(zip(eigenvalues, dets), key=lambda a:-abs(a[0])):
		print "%.6f %s" % (ev, det)

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print >> sys.stderr, """\
			Usage: %s [configurations file] [coefficients file]

				[configurations file]: Gaussian .log-file containing a list of all configurations.
				                       In case of large CASSCF-calculations, this requirs
				                           iop(4/43=3)
				                       in the input.
				[coefficients file]: Gaussian .log-file containing all eigenvalues in the format
				                     given when using iop(5/33=2).
		""" % sys.argv[0]
		sys.exit(1)

	main(sys.argv[1], sys.argv[2])
