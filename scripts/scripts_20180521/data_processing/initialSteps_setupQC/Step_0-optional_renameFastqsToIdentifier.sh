# Want to give the Baja samples consistent identifiers

#AB1GN1SS01 --> 168_Elut_BAJ_TS2
#AB1GN1SS01  --> 169_Elut_BAJ_R1
#AB1GN1SS03 --> 170_Cfam_CONTROL_LANGEDOG

# you can use the *magical* rename command to do this!
# usage:
# rename [pattern] [replacement pattern] [files]

rename AB1GN1SS01 168_Elut_BAJ_TS2 AB1GN1SS01*

rename AB1GN1SS02 169_Elut_BAJ_R1 AB1GN1SS02*

rename AB1GN1SS03 170_Cfam_CONTROL_LANGEDOG AB1GN1SS03*

