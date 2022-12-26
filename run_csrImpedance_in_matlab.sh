#!/bin/bash

args=("$@")

csrImpedance Z_CSR_pp.sdds -height=${args[0]} -radius=${args[1]} -frequencyLimit=maximum=${args[2]} -n=10

sdds2plaindata Z_CSR_pp.sdds Z_CSR_pp.txt -outputMode=ascii -noRowCount "-separator= " -col=k -col=ZReal -col=ZImag

exit 0