# catie bullen
# dec 2022
# quick fix for DEV-1578 until long-term solution implemented
# usage: python dcf_manifest_acl_fix.py <file>

import json
import pandas as pd
import sys


# set up error log file and error stream
err_out = "./acl_fix_error_log.txt"
sys.stderr = open(err_out, "w")

# list of offender projects
# note these are only released projects in GDC, more projects exist that are not yet released to public
offenders = [
    "CGCI-BLGSP",
    "TARGET-CCSK",
    "CTSP-DLBCL1",
    "CGCI-HTMCP-CC",
    "TARGET-ALL-P2",
    "TARGET-ALL-P1",
    "TARGET-AML",
    "TARGET-WT",
    "TARGET-OS",
    "TARGET-RT",
    "TARGET-NBL",
]

# project-level ACLs (want to REMOVE these from the manifest list object under acl column)
# this is for reference

project_dict = {
    "CGCI-BLGSP": "phs000527",
    "TARGET-CCSK": "phs000466",
    "CTSP-DLBCL1": "phs001184",
    "CGCI-HTMCP-CC": "phs000528",
    "TARGET-ALL-P2": "phs000464",
    "TARGET-ALL-P1": "phs000463",
    "TARGET-AML": "phs000465",
    "TARGET-WT": "phs000471",
    "TARGET-OS": "phs000468",
    "TARGET-RT": "phs000470",
    "TARGET-NBL": "phs000467",
}

# program-level ACLs (want to KEEP/USE these in the manifest list object under acl column)

program_dict = {
    "CGCI-BLGSP": "phs000235",
    "TARGET-CCSK": "phs000218",
    "CTSP-DLBCL1": "phs001175",
    "CGCI-HTMCP-CC": "phs000235",
    "TARGET-ALL-P2": "phs000218",
    "TARGET-ALL-P1": "phs000218",
    "TARGET-AML": "phs000218",
    "TARGET-WT": "phs000218",
    "TARGET-OS": "phs000218",
    "TARGET-RT": "phs000218",
    "TARGET-NBL": "phs000218",
}

# function to apply to dataframe for correct acl


def switcheroo(project_id):
    return [program_dict[project_id]]


# read in data frame

df = pd.read_csv(sys.argv[1], sep="\t")

# separate records by whether or not belong to problem projects or not

non_offender_df = df[~df.project_id.isin(offenders)]

offender_df = df[df.project_id.isin(offenders)]

# apply function to problem records

offender_df["acl"] = offender_df["project_id"].apply(switcheroo)

# rejoin  data sets

df_v2 = pd.concat([non_offender_df, offender_df])

# write back to file

df_v2.to_csv("ACL_update_" + sys.argv[1], sep="\t", index=False)
