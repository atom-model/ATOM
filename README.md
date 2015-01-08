# ATOM
* * *

###  What is this?
ATOM (Atmospheric and Ocean Model) is a fast climate model.
ADD MORE DETAILS

##### Requirements:
- Global bathymetry/topography grid in 1° x 1° spacing -- paleotopography/bathymetry grids (Smith et al. 1994; Golonka et al. 1997) between 140 - 0 Ma are included here, created using agegrid rev.210 (or Earthbyte 2013.2.rot)
- Present day surface temperature: included based on NASA
- Present day precipitation: included from NASA
- Present day salinity: included from NASA

This bitbucket repository requires mercurial ([http://mercurial.selenic.com/downloads]
(http://mercurial.selenic.com/downloads)). Macports:  sudo port install mercurial

* * *

### How do I download ATOM?
In terminal: `hg clone https://yourusername@bitbucket.org/nickywright/atom`


### How do I add or change stuff?

1. In terminal, navigate to the atom directory
2.
    - To modify a pre-exisiting file: `hg commit -m 'YOUR COMMENT HERE'`. You must add a comment to commit
    - To add a new file:
          - `hg st` to list new files - files that haven't been added to the repo will have a preceding `?`
          - `hg add FILENAME` (to add an individual file) or `hg add` (to add *all* files including invisible files)
          - (Optional): hg st - files that you want to commit will now have a preceding `A`
          - `hg commit -m 'YOUR COMMENT HERE' `. You must add a comment to commit
3. Push edits to bitbucket by typing in `hg push` - this will ask for your bitbucket password

For more information, please look here: [http://mercurial.selenic.com/guide](http://mercurial.selenic.com/guide)

### How do I update my local version ###
1. In terminal: `hg pull` (This doesn't modify your files though!)
2. `hg update`  (or maybe hg update --check)- this updates your files, based on the changes it downloaded in `hg pull` (I think)

### How do I download a previous commit?
In your browser window (Safari, Firefox, Chrome), type:

           https://**YOUR USERNAME**@bitbucket.org/nickywright/atom/get/**COMMITNUMBER**.tar.gz

COMMITNUMBER can be found under the Commit heading on https://bitbucket.org/nickywright/atom/commits/all

* * *

### Summary of changes
07 Jan 2015: Version 1.0
