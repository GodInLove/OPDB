def extract(infile, outfile):
    f = open(infile, 'r')
    out_f = open(outfile, 'w')
    out_f.write("fwd,rev\n")
    while True:
        line = f.readline().strip()
        if not line:
            break
        li = line.split()
        out_f.write(li[2] + "," + li[3] + "\n")

    f.close()
    out_f.close()