import sys
import math

# run with::
# ./bin/vg view -a alt.mega.gam | grep -o "\"length\": [-]*[0-9]*" | cut -f 2 -d ":" | python scripts/calc_insert.py

def calc_mean(vals):
    return float(sum(vals)) / float(len(vals))

def calc_sd(vals, mean):
    sq_err_sum = 0.0
    for i in vals:
        sq_err_sum += (float(i) - mean) ** 2

    return  math.sqrt( (1.0 / float(len(vals))) * sq_err_sum)

if __name__ == "__main__":

    vals = []
    for line in sys.stdin:
        x = float(line.strip())
        if x > 0:
            vals.append(x);
        else:
            vals.append(-1.0 * x)
        if len(vals) == 1000:
            mean = calc_mean(vals)
            sd = calc_sd(vals, mean)
            print("Mean: ", mean, ",  SD: ", sd)
            del vals[:]
