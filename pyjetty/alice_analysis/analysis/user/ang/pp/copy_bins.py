''' Program to copy theory bins onto another set of predictions
    Ezra Lesser (elesser@berkeley.edu), Spring 2021
'''

import os       # File operations
from pathlib import Path
import bisect   # Some sorted list operations


# User-defined parameters
old_path = "/home/ezra/Archive/"     # Old data points to merge into
new_path = "/home/ezra/firstbin/"    # New data points to insert into old
create_path = "/home/ezra/merged/"   # Path for newly created data

pT_list = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, \
           100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
jetR_list = [0.2, 0.4]
beta_list = [1.5, 2, 3]


# Define helper functions
def get_dir(base, gr, jetR, beta, pTmin, pTmax,
            check_exists=False, create_new=False):

    f_dir = os.path.join(
        base, gr + "ALICE_R%s" % str(jetR).replace('.', ''),
        "beta%s" % str(beta).replace('.', 'p'),
        "pT%i_%i" % (pTmin, pTmax))
    if check_exists and not os.path.exists(f_dir):
        raise ValueError("Directory does not exist: %s" % f_dir)
    if create_new:
        Path(f_dir).mkdir(parents=True, exist_ok=True)
    return f_dir


def get_val_li(file):
    lines = [line for line in file.read().split('\n') if line]
    x_li = [float(line.split()[0]) for line in lines]
    y_li = [float(line.split()[1]) for line in lines]
    return (x_li, y_li)


def create_merged_files(new_dir, old_dir, create_dir):
    # Loop through all files in the directory
    for l in [0, 1, 2]:
        for m in [0, 1, 2]:
            for n in [0, 1, 2]:
                # Skip non-existant scale variations
                if (0 in (l, m, n) and 2 in (l, m, n)):
                    continue
                create_merged_file_lmn(
                    new_dir, old_dir, create_dir, l, m, n)


def create_merged_file_lmn(new_dir, old_dir, create_dir, l, m, n):
    # Initialize & load data points from both files
    x_li_new = None; y_li_new = None;
    x_li_old = None; y_li_old = None;
    with open(os.path.join(new_dir, "%i%i%i.dat" % (l, m, n)), 'r') as f:
        x_li_new, y_li_new = get_val_li(f)
    with open(os.path.join(old_dir, "%i%i%i.dat" % (l, m, n)), 'r') as f:
        x_li_old, y_li_old = get_val_li(f)

    # Find correct place to do the splice
    min_i = bisect.bisect_left(x_li_old, x_li_new[0])
    max_i = bisect.bisect(x_li_old, x_li_new[-1], lo=min_i)

    # Renormalize the new points to match the expected value of
    #     the last old point that was removed
    if min_i != max_i:  # simultaneously requires max_i > 0
        x_old = x_li_old[max_i-1]
        y_old = y_li_old[max_i-1]
        try:  # if x_old in x_li_new
            # Renormalize purely by one point
            y_new = y_li_new[x_li_new.index(x_old)]
            y_li_new = [y * y_old / y_new for y in y_li_new]
        except ValueError:  # x_old not in x_li_new
            # Renormalize by avg of two points
            min_i_new = bisect.bisect_left(x_li_new, x_old)
            avg = (y_li_new[min_i_new-1] + y_li_new[min_i_new]) / 2.
            y_li_new = [y * y_old / avg for y in y_li_new]

    x_li = x_li_old[0:min_i] + x_li_new + x_li_old[max_i:-1]
    y_li = y_li_old[0:min_i] + y_li_new + y_li_old[max_i:-1]

    with open(os.path.join(create_dir, "%i%i%i.dat" % (l, m, n)), 'w') as f:
        f.writelines(["{:.4e}  {:.4e}\n".format(x_li[i], y_li[i]) for i in range(len(x_li))])


# Loop through all directories
for jetR in jetR_list:
    print("Working on R =", jetR, "...")
    for i, pTmin in enumerate(pT_list[:-1]):
        pTmax = pT_list[i+1]
        for beta in beta_list:
            for gr in ["", "gr_"]:

                # Load file directories
                new_dir = get_dir(
                    new_path, gr, jetR, beta, pTmin, pTmax, check_exists=True)

                old_dir = get_dir(
                    old_path, gr, jetR, beta, pTmin, pTmax, check_exists=True)

                create_dir = get_dir(
                    create_path, gr, jetR, beta, pTmin, pTmax, create_new=True)

                create_merged_files(new_dir, old_dir, create_dir)
