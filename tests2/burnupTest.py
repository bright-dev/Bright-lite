import re
import os
import subprocess

main_dir = os.getcwd()
subprocess.call("rm -r outBUD", shell = True)
enrich = 20.
while (enrich < 45):
    os.mkdir(str(enrich)+"%")
    subprocess.call("cp Bright-lite " + str(enrich) + "%", shell = True)
    f = open("inputFile.temp", 'r')
    lines = f.readlines()
    f.close()
    lines[0] = re.sub('ENRICHMENT', str(enrich), lines[0])
    lines[1] = re.sub('ENRICHMENT2', str(1000. - enrich), lines[1])
    g = open("inputFile.txt", 'w')
    for line in lines:
        g.write(line)
    g.close()
    os.chdir(str(enrich)+"%")
    subprocess.call("./Bright-lite", shell= True)
    os.chdir(main_dir)
    enrich = enrich + 0.5

   
     
