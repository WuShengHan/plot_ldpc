#!/usr/bin/python

import sys, string, re

name = 'ber2snr.py'
fin = open(name, "r")
fout = open("gen_%s.py" % name, "w")

fout.write("#!/usr/bin/python\n\nimport sys, string, re\n\n")
fout.write("fout = open(\"tmp_%s\", \"w\")\n" % name)
fout.write("\n")

for line in iter(fin):
    inline = line.splitlines()
    oline = inline[0]
    oline = oline.replace('\\','\\\\')
    oline = re.sub("\"","\\\"", oline)
    #print oline
    oline = "fout.write(\"%s\\n\")\n" % oline
    fout.write(oline)

fout.write("fout.close()\n")
fin.close()
fout.close()


#line = ""
#line = line + """
#always@(fail_num)
#begin
#   $display($time, ", Error: num compare mismatch!");
#   #1000;  $finish;
#end
#
#"""
#print >>fp_log, """
#always@(fail_loc)
#begin
#   $display($time, ", Error: loc compare mismatch!");
#   #1000;  $finish;
#end
#"""
#fp_log.write(line)

