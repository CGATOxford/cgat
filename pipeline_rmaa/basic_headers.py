
from subprocess import Popen, PIPE
import os
from time import sleep
import commands

def run_cmd(my_cmd):
   p=Popen(my_cmd, shell=True, stdin = PIPE, stderr = PIPE)
   p.stdin.close()
   output=p.stderr.read()
   sts = os.waitpid(p.pid, 0)
   if 'error' in output:
      print "could not contact SGE, retrying resubmission in 5 seconds"
      print my_cmd
      sleep(5)
      run_cmd(my_cmd)
   return

def run_cmd_output(my_cmd):
   p=Popen(my_cmd, shell=True, stdin = PIPE, stdout = PIPE, stderr = PIPE)
   p.stdin.close()
   output=p.stdout.read()
   eout=p.stderr.read()
   sts = os.waitpid(p.pid, 0)
   if 'error' in eout:
      print "could not contact SGE, retrying resubmission in 5 seconds"
      print my_cmd
      sleep(5)
      output=run_cmd_output(my_cmd)
   return output

def run_cmd_qsub(my_cmd):
   while True:
      sleep(5)
      if 'has been submitted' in commands.getoutput(my_cmd):
         break
   return

