# usu-biol4750
A place for students to go to look at records from live coding examples and notes.

# Logging into the CHPC server:  
One way to get to your CHPC user interface and have access to your home directory and group storage directory:
https://ondemand.chpc.utah.edu/pun/sys/dashboard/files/fs/uufs/chpc.utah.edu/common/home/u6036559  
...replacing uXXXXXXX with your U of U user ID.  

# Connecting to RStudio server on CHPC:
https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sys/rstudio_server_app/session_contexts/new        
**R version:** R 4.4.0 Geospatial packages   
**Cluster:** notchpeak      
**Account:Partition** saarman-np:saarman-shared-np (very important!!! this allows multiple jobs to run simultaneously, if you don't use this partition you will block entry for everyone else!)      
**Number of cores:** 1-4 (there are 32 total)   
**Number of hours:** 100 (336 is the max, definitely put more than 24)   
**Memory per job:** 128G (1000 GB total is the max across all of the 32 cores + all users, but good to stay well below half of the limit)  

**You can access a running session by going to https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sessions 

# Helpful links:  
  
- To start a new R session, use this link: <https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sys/rstudio_server_app/session_contexts/new>
     
    R version: **R 4.4.0 Geospatial packages** #changed from R 4.0.3 !    
    Cluster: notchpeak  
    Account and partition: usu-biol4750:notchpeak-shared-freecycle   
    Number of cores (per node): 4   
    Number of hours: 72  
    Memory per job in GB: 16
   
- To reconnect to an existing session, you can go to <https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sessions>

- Lab Manual: <https://bookdown.org/hhwagner1/LandGenCourse_book/>

- Poll everywhere: <https://PollEv.com/norahsaarman260>

- Textbook: <https://libcat.lib.usu.edu/record=b3469639~S11>  
  
# Example Lab Portfolios  
- This is Where Norah will save notes from live coding examples in class, and serves as an example of the portfolio entry you will turn in for all/nothing credit in Canvas.
 - "Save As..." to create a copy of this file with your name/initials, and then fill in the gaps to keep track of your coding, thoughts, and learning.
   
## Lab 1: Review of R Skills 9/3/2024
 - lab01-NPS-example.Rmd (lab01-NPS-example.pdf is the knitted version)
 - Following along from <https://bookdown.org/hhwagner1/LandGenCourse_book/basic-r.html>
 - Create an Rmd log entry that provides the code from:
 - **a. Section 2.1 Basic R Programming** (you can mostly copy and paste and add comments)
 - **b. Section 2.2 R Graphics** (you can mostly copy and paste and add comments)
 - At the end, export your work as .Rmd/.txt/.md/etc from RStudio, and submit to Canvas for credit and feedback. Ask for help if you need it!

## Lab 2: Importing Genetic Data 9/10/2023
 - lab02-NPS-example.Rmd  (lab02-NPS-example.pdf is the knitted versions)  
 - Following along from https://bookdown.org/hhwagner1/LandGenCourse_book/Week1.html
 - Create an Rmd log entry that provides the code from:
 - **a. Section 4.3 Worked Example** (you can mostly copy and paste and add comments, 4.3.5 and on is OPTIONAL)
 - **b. Section 4.4 R exercise Week 1** (you will need to write your own code using knowledge from Section 4.3)
 - At the end, **knit it to a pdf or md file**, export from RStudio, and submit to Canvas for credit and feedback. If there are error messages during knitting, check the ‘R Markdown’ tab for the code line number and try to fix it. Ask for help if you need it!

## Lab 3: Genetic Diversity 9/17/2024  
 - lab03-NPS-example.Rmd (lab03-NPS-example.pdf is the knitted versions)
 - Following along from https://bookdown.org/hhwagner1/LandGenCourse_book/Week3.html
 - Create an Rmd log entry that provides the code from:
 - **a. Section 6.3 Worked Example** (you can mostly copy and paste and add comments)
 - **b. Section 6.4 R exercise Week 3** (you will need to write your own code using knowledge from Section 6.3)
 - At the end, **knit it to a pdf or md file**, export from RStudio, and submit to Canvas for credit and feedback. If there are error messages during knitting, check the ‘R Markdown’ tab for the code line number and try to fix it. Ask for help if you need it!

## Lab 4: Spatial Data 10/3/2024 
 - lab04-NPS-example.Rmd (lab04-NPS-example.pdf is the knitted versions)
 - Following along from https://bookdown.org/hhwagner1/LandGenCourse_book/Week2.html
 - Create an Rmd log entry that provides the code from:
 - **Section 5.4 R exercise Week 2** (you will need review applicable code from Section 5.3 to write your own code)
 - Optional: Add to this log entry with 5.5 Spatial Data – Bonus Vignette: 'sf' package, plotting categorical maps
 - At the end, **knit it to a pdf or md file**, export from RStudio, and submit to Canvas for credit and feedback. If there are error messages during knitting, check the ‘R Markdown’ tab for the code line number and try to fix it. Ask for help if you need it!
   
## Lab 5: Spatial Statistics 10/3/2024 
 - lab05-NPS-example.Rmd (lab05-NPS-example.pdf is the knitted versions)
 - Following along from https://bookdown.org/hhwagner1/LandGenCourse_book/Week5.html
 - Create an Rmd log entry that provides the code from:
 - **Section 8.5 R exercise Week 5** (you will need review applicable code from Section 5.3 to write your own code)
 - At the end, **knit it to a pdf or md file**, export from RStudio, and submit to Canvas for credit and feedback. If there are error messages during knitting, check the ‘R Markdown’ tab for the code line number and try to fix it. Ask for help if you need it!
 - **Optional Lab 5 bonus: Moran's I** = Section 8.3 R worked example last section - Moran's I (this code uses the package 'spdep')

## Lab 6: Resistance Surfaces 11/11/2024
 - lab06-NPS-example.Rmd (lab06-NPS-example.pdf is the knitted versions)
 - Following along from https://bookdown.org/hhwagner1/LandGenCourse_book/Week13.html
 - Create an Rmd log entry that provides the code from:
 - **Section 13.1 R worked example**
 - At the end, **knit it to a pdf or md file**, export from RStudio, and submit to Canvas for credit and feedback. If there are error messages during knitting, check the ‘R Markdown’ tab for the code line number and try to fix it. Ask for help if you need it!

