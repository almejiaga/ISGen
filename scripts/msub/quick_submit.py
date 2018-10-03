import sys, os
import subprocess

print sys.argv

header = os.path.expanduser(sys.argv[1])
outdir = os.path.expanduser(sys.argv[2])
scriptdir = os.path.expanduser(sys.argv[3])
command = ' '.join(sys.argv[4:])

job_file = os.path.join(outdir, 'job_script.sh')
subprocess.call(['cp', header, job_file])

with open(job_file, 'a') as f:
    f.write('#PBS -o ' + outdir + '\n')
    f.write('#PBS -e ' + outdir + '\n')
    f.write('\n')
    f.write('cd ' + scriptdir + '\n')
    f.write(command + '\n')

subprocess.call(['qsub', job_file])
