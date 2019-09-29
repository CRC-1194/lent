#!/usr/bin/env python3

import testReportCore as trc

def main():
    hardConditions = ["consistent_normals", "unit_length"]
    trc.compile_hard_condition_violation_tex_table(hardConditions)

if __name__ == "__main__":
    main()
