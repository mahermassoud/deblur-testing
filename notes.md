# Sample intersection!!
In EMP release 1. Everything in 150nt is in 100nt and 90nt except for these 2 samples
                                                 sample        count
34                                         1064.G.CV298        5.000
2229.W2.N13.EH1.Thomas.CMB.Seaweed.lane5.NoIndex.L005

# 4/17/18 TODO
1. Subsample EMP and save 10 - 100 samples
2. Set up do you can plot runtime vs # samples (do in pbs script?)
    - keep track of #features
3. Run and plot runtime v #samples
4. Maybe do same for features?

# High priority
1. Make dotfile for bash config, use miniconda
2. Put it on cluster
3. Make interactive job
4. simple qsub / pbs script
5. how to request 3 processors and 24gb of memory, walltime of 2 hr (kill after 2 hr)
    - know how to kill a job thats running or in the queue
6. how to automaticaly setup a temp directory on panfs
    - when done move it to home, delete
7. Run unit tests in command line    
8. Write script to see if num_cores param makes a difference!



# Low priority
- do eg.s in python console format
- try with normalized data, non euclidian, divide vector by sum of vector
- get rid of percents business
- Fix Demux test in methods by using small samples
- look @ calculating trim length
- input validity checks
- make clps_df flow better
- lint
- make -o default to . instead of requiring
- get rid of docstrings from click commands and improve interface, --help
- more metrics!! eg. seq depth
     -Look @ distribution of difference-per-featureplot taxa that got dropped out
     -look @ distribution of top collapsed features
     -if it is non-uniform, that would suggest we are losing info
        - ALSO need to look at weight.
        - eg. if 90% of an OTU's counts are collapsed is that bad?
        - plot percent of counts in otu that are collapsed?
     -Straight up jac/bc/mantel between biom tables
     -sequencing depth

# Can't do yet
- test against different environments eg. fecal, skin, soil

# Ideas
- append sample name in pre_post

# Questions
- When they said plot taxa, did they mean the actual species?

# Meeting notes 4/19/18
- make a dotfile github bash_profile, copy dans
- depenencies, download miniconda, make conda inv
    1) make environment
    2) make interactive job, understand when in and out of the job
    3) qsub/pbs script, something simple like submitting hostname
        - dont worry about amount of hardware
    4) how to request 3 processors and 24gb of memory, walltime of 2 hr (kill after 2 hr)
        - know how to kill a job thats running or in the queue
    5) how to automaticaly setup a temp directory on panfs
        - when done move it to home, delete
        
- google torque examples
- screen / tmux for running mutliple different shells
- run unit tests on new environment-- figure out how to do on terminal

# Notes 5/6/18
modify ranalysis so you can run it on different trim lengths