# R3B_Neutron_Generator

Use this file to generate input file to be used with the R3BRoot simulation. Allows to generate 1n/2n decay events as well as different decay mechanisms in case of the 2n decay.

Different options can be selected by modifying the file and uncommenting/commenting some parts.

Download :

>git clone https://github.com/aldros/R3B_Neutron_Generator.git

>cd R3B_Neutron_Generator

>git checkout dev

To run : 

>root

>.L NeutronDecayGenerator.cxx++

>EventGenerator_Ndecay("fileName",10000) //arguments are fileName and number of events to consider
