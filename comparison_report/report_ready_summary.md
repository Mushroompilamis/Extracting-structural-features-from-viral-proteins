# Report-ready capsid contact comparison

## Analysis design

ChimeraX `crystalcontacts` outputs were compared for the WT/control and mutant capsid models. Residue numbers were normalized so 1BMV-style observed numbers such as 3159 are interpreted as sequence residue 159. Thus, the report can refer to T159A, while the ChimeraX command remains `swapaa /2:3159 ala`.

## Critical quality check

The WT and mutant contact files must be generated using the same ChimeraX protocol: same structure state, same assembly/ASU choice, same distance cutoff, and same `crystalcontacts` settings. If WT was calculated from a different model state than the mutants, very large global differences may reflect protocol differences rather than mutation effects.

## Side-by-side global metrics

```text
sample  raw_contacts  unique_observed_contacts  unique_sequence_contacts  inter_chain_or_copy_contacts  symmetry_related_contacts  hydrogen_bond_like_contacts  hydrophobic_contacts  electrostatic_contacts  aromatic_contacts  mean_distance  median_distance  contacts_le_3_5A  contacts_le_4A  total_buried_sasa  mean_buried_sasa  max_buried_sasa  buried_residue_entries
    WT          2227                      2053                      1084                          2227                        451                          322                   242                     112                 52          3.817             3.81               724            1354            940.331             1.353             20.8                     695
 T159A           749                       684                       351                           749                        131                           99                    82                      38                 12          3.847             3.90               227             416            415.468             1.468             14.2                     283
 T155A           753                       686                       352                           753                        135                          111                    80                      38                 12          3.833             3.88               231             422            464.588             1.636             19.0                     284
```

## Key changes vs WT/control

```text
sample                       metric  control_value  sample_value     delta  percent_change
 T159A     unique_sequence_contacts       1084.000       351.000  -733.000          -67.62
 T159A inter_chain_or_copy_contacts       2227.000       749.000 -1478.000          -66.37
 T159A    symmetry_related_contacts        451.000       131.000  -320.000          -70.95
 T159A  hydrogen_bond_like_contacts        322.000        99.000  -223.000          -69.25
 T159A         hydrophobic_contacts        242.000        82.000  -160.000          -66.12
 T159A            total_buried_sasa        940.331       415.468  -524.863          -55.82
 T155A     unique_sequence_contacts       1084.000       352.000  -732.000          -67.53
 T155A inter_chain_or_copy_contacts       2227.000       753.000 -1474.000          -66.19
 T155A    symmetry_related_contacts        451.000       135.000  -316.000          -70.07
 T155A  hydrogen_bond_like_contacts        322.000       111.000  -211.000          -65.53
 T155A         hydrophobic_contacts        242.000        80.000  -162.000          -66.94
 T155A            total_buried_sasa        940.331       464.588  -475.743          -50.59
```

## Interpretation

### T159A vs WT

- unique sequence contacts changed from 1084.0 to 351.0 (-733.000, -67.62%).
- inter chain or copy contacts changed from 2227.0 to 749.0 (-1478.000, -66.37%).
- symmetry related contacts changed from 451.0 to 131.0 (-320.000, -70.95%).
- hydrogen bond like contacts changed from 322.0 to 99.0 (-223.000, -69.25%).
- hydrophobic contacts changed from 242.0 to 82.0 (-160.000, -66.12%).
- total buried sasa changed from 940.331 to 415.468 (-524.863, -55.82%).

This suggests altered interface/contact organization in the mutant model. The result should be interpreted as an in silico structural-contact comparison, not direct experimental proof that capsid assembly fails.

### T155A vs WT

- unique sequence contacts changed from 1084.0 to 352.0 (-732.000, -67.53%).
- inter chain or copy contacts changed from 2227.0 to 753.0 (-1474.000, -66.19%).
- symmetry related contacts changed from 451.0 to 135.0 (-316.000, -70.07%).
- hydrogen bond like contacts changed from 322.0 to 111.0 (-211.000, -65.53%).
- hydrophobic contacts changed from 242.0 to 80.0 (-162.000, -66.94%).
- total buried sasa changed from 940.331 to 464.588 (-475.743, -50.59%).

This suggests altered interface/contact organization in the mutant model. The result should be interpreted as an in silico structural-contact comparison, not direct experimental proof that capsid assembly fails.

## Target residue neighborhood

```text
sample target_sequence_residue target_observed partner_sequence_residue partner_observed  distance   interaction_type  inter_chain_or_copy  symmetry_related
 T155A                  THR159 THR3159_2_ncs14                    ASN35  ASN3035_2_ncs15      1.42 hydrogen_bond_like                 True              True
 T155A                  THR159 THR3159_2_ncs14                    ASN35  ASN3035_2_ncs15      1.42 hydrogen_bond_like                 True              True
 T155A                  THR159  THR3159_2_ncs8                   THR159   THR3159_2_ncs4      2.26 hydrogen_bond_like                 True              True
 T155A                  THR159  THR3159_2_ncs4                   THR159   THR3159_2_ncs8      2.26 hydrogen_bond_like                 True              True
 T155A                  THR159 THR3159_2_ncs14                    LYS34  LYS3034_2_ncs15      2.74 hydrogen_bond_like                 True              True
 T155A                  THR159 THR3159_2_ncs14                    LYS34  LYS3034_2_ncs15      2.74 hydrogen_bond_like                 True              True
 T155A                  THR159 THR3159_2_ncs15                    ASN35  ASN3035_2_ncs14      2.76 hydrogen_bond_like                 True              True
 T155A                  THR159 THR3159_2_ncs15                    ASN35  ASN3035_2_ncs14      2.76 hydrogen_bond_like                 True              True
 T155A                  THR159  THR3159_2_ncs4                    LYS34   LYS3034_2_ncs8      3.01 hydrogen_bond_like                 True              True
 T155A                  THR159  THR3159_2_ncs4                    LYS34   LYS3034_2_ncs8      3.01 hydrogen_bond_like                 True              True
 T155A                  THR159  THR3159_2_ncs4                    ASN35   ASN3035_2_ncs8      3.15 hydrogen_bond_like                 True              True
 T155A                  THR159  THR3159_2_ncs4                    ASN35   ASN3035_2_ncs8      3.15 hydrogen_bond_like                 True              True
 T155A                  THR159 THR3159_2_ncs15                    LYS34  LYS3034_2_ncs14      3.55              other                 True              True
 T155A                  THR159 THR3159_2_ncs15                    LYS34  LYS3034_2_ncs14      3.55              other                 True              True
 T155A                  THR159 THR3159_2_ncs15                   THR159  THR3159_2_ncs14      3.63              other                 True              True
 T155A                  THR159 THR3159_2_ncs14                   THR159  THR3159_2_ncs15      3.63              other                 True              True
 T155A                  THR159  THR3159_2_ncs8                    ASN35   ASN3035_2_ncs4      3.76              other                 True              True
 T155A                  THR159  THR3159_2_ncs8                    ASN35   ASN3035_2_ncs4      3.76              other                 True              True
 T155A                  THR159  THR3159_2_ncs4                   PRO160   PRO3160_2_ncs8      4.33              other                 True              True
 T155A                  THR159  THR3159_2_ncs4                   PRO160   PRO3160_2_ncs8      4.33              other                 True              True
 T155A                  THR159  THR3159_2_ncs8                    LYS34   LYS3034_2_ncs4      4.35              other                 True              True
 T155A                  THR159  THR3159_2_ncs8                    LYS34   LYS3034_2_ncs4      4.35              other                 True              True
    WT                  THR159 THR3159_2_ncs14                    ASN35  ASN3035_2_ncs15      1.42 hydrogen_bond_like                 True              True
    WT                  THR159 THR3159_2_ncs14                    ASN35  ASN3035_2_ncs15      1.42 hydrogen_bond_like                 True              True
    WT                  THR159  THR3159_2_ncs8                   THR159   THR3159_2_ncs4      2.26 hydrogen_bond_like                 True              True
    WT                  THR159  THR3159_2_ncs4                   THR159   THR3159_2_ncs8      2.26 hydrogen_bond_like                 True              True
    WT                  THR159 THR3159_2_ncs14                    LYS34  LYS3034_2_ncs15      2.74 hydrogen_bond_like                 True              True
    WT                  THR159 THR3159_2_ncs14                    LYS34  LYS3034_2_ncs15      2.74 hydrogen_bond_like                 True              True
    WT                  THR159 THR3159_2_ncs15                    ASN35  ASN3035_2_ncs14      2.76 hydrogen_bond_like                 True              True
    WT                  THR159 THR3159_2_ncs15                    ASN35  ASN3035_2_ncs14      2.76 hydrogen_bond_like                 True              True
    WT                  THR159  THR3159_2_ncs4                    LYS34   LYS3034_2_ncs8      3.01 hydrogen_bond_like                 True              True
    WT                  THR159  THR3159_2_ncs4                    LYS34   LYS3034_2_ncs8      3.01 hydrogen_bond_like                 True              True
    WT                  THR159  THR3159_2_ncs4                    ASN35   ASN3035_2_ncs8      3.15 hydrogen_bond_like                 True              True
    WT                  THR159  THR3159_2_ncs4                    ASN35   ASN3035_2_ncs8      3.15 hydrogen_bond_like                 True              True
    WT                  THR159 THR3159_2_ncs15                    LYS34  LYS3034_2_ncs14      3.55              other                 True              True
    WT                  THR159 THR3159_2_ncs15                    LYS34  LYS3034_2_ncs14      3.55              other                 True              True
    WT                  THR159 THR3159_2_ncs15                   THR159  THR3159_2_ncs14      3.63              other                 True              True
    WT                  THR159 THR3159_2_ncs14                   THR159  THR3159_2_ncs15      3.63              other                 True              True
    WT                  THR159  THR3159_2_ncs8                    ASN35   ASN3035_2_ncs4      3.76              other                 True              True
    WT                  THR159  THR3159_2_ncs8                    ASN35   ASN3035_2_ncs4      3.76              other                 True              True
    WT                  THR159  THR3159_2_ncs4                   PRO160   PRO3160_2_ncs8      4.33              other                 True              True
    WT                  THR159  THR3159_2_ncs4                   PRO160   PRO3160_2_ncs8      4.33              other                 True              True
    WT                  THR159  THR3159_2_ncs8                    LYS34   LYS3034_2_ncs4      4.35              other                 True              True
    WT                  THR159  THR3159_2_ncs8                    LYS34   LYS3034_2_ncs4      4.35              other                 True              True
```

## Lost/gained sequence-level contacts

### T159A

- Lost contacts vs WT: 736
- Gained contacts vs WT: 3

Examples of lost contacts:
- ALA116 -- ALA22
- ALA116 -- ASP23
- ALA116 -- PRO24
- ALA116 -- VAL25
- ALA117 -- ARG151
- ALA117 -- ASN114
- ALA119 -- VAL112
- ALA121 -- PRO32
- ALA140 -- ARG60
- ALA140 -- LYS63
- ALA140 -- PHE61
- ALA140 -- PHE62
- ALA140 -- THR22
- ALA145 -- THR23
- ALA167 -- VAL206

### T155A

- Lost contacts vs WT: 732
- Gained contacts vs WT: 0

Examples of lost contacts:
- ALA116 -- ALA22
- ALA116 -- ASP23
- ALA116 -- PRO24
- ALA116 -- VAL25
- ALA117 -- ARG151
- ALA117 -- ASN114
- ALA119 -- VAL112
- ALA121 -- PRO32
- ALA140 -- ARG60
- ALA140 -- LYS63
- ALA140 -- PHE61
- ALA140 -- PHE62
- ALA140 -- THR22
- ALA145 -- THR23
- ALA167 -- VAL206

## Suggested sentence for the report

The in silico T→A mutations reduced several interface descriptors, including inter-chain/copy contacts, symmetry-related contacts, hydrogen-bond-like contacts, and buried surface area compared with the WT control. This supports the hypothesis that the selected threonine residues contribute to inter-subunit stabilization in the 1BMV capsid, although experimental validation would be required to confirm an assembly defect.
