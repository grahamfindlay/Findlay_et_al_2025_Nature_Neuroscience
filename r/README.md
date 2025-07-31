# R evironments and IDEs on `tononi-2`. 
Using system `R` on `tononi-2` is terrible! 
The combination of the system-wide site library (that you cannot write to without root), and a personal library, never "just works". Don't bother! 
Use `pixi` or `mamba` (see below). Unfortunately [`rig`](https://github.com/r-lib/rig) requires admin. 

## Using `pixi` to create a user R on `tononi-2`
Getting VS Code's R extension to play nicely with R virtual environments is hell. 
Using `pixi`'s `global` installation feature is a compromise. Better than using system R. 
A good [source](https://wanggroup.org/productivity_tips/vscode-setup.html). 
```{bash}
pixi global install r-base
```
This creates an environment also called `r-base` at `~/.pixi/envs/r-base`. 
I don't love that this is implicit. But it works, and `.libPaths()` doesn't find system R. 
```{bash}
pixi global install --environment r-base r-arrow r-dplyr etc.
```
This works too, but system `radian` is usually fine. 
```{bash}
pixi global install --environment r-base radian
```
Make sure to log out and back in after installing `r-base` and `radian`, or else the latter will continue to point to system R. 

Unfortuantely, there is no `conda-forge` feedstock for [`vscDebugger`](https://github.com/ManuelHentschel/vscDebugger). 

Pros:
- Because many R ecosystem tools actually *expect* you to be using a system/global R installation, rather than a virtual environment (/sigh), this generally works well and avoids a lot of headaches. 
Cons: 
- Unforunately, at time of writing, Positron does not currently support discovery of `pixi` R environments, even though it seems like they should work according to the Positron docs. You can give Positron the path to your `pixi` R, and it still will refuse to use it. 


## Using `mamba` to create a user R on `tononi-2`. 
Simple. Use [`miniforge`](https://github.com/conda-forge/miniforge), and the included `environment.yml`, to install `R`:
```{bash}
conda env create -f environment.yml
```

Pros: 
- Works with Positron.
- Provides an actual virtual environment. 
Cons: 
- See above note about many R tools actually expecting system R.
- There is no other case in which `mamba` > `pixi`. 

# Using VS Code
For code intelligence to work with your functions, you MUST open the package directory (`f25a/`) as a workspace!

*If* you are working from a parent directory, be mindful that VS Code's R extension will set your initial working directory differently than RStudio or Positron. You may want to stick a `.Rprofile` at the root of your VS Code project.
For example:
```{r}
setwd("/path/to/f25a")
```
This will ensure that this is your initial working directory, and that `here` beings moving up the filesystem hierarchy from there. 

# Using Positron
Getting Positron to recognize any non-system R seems to be challenging. 
Positron does not recognize even a globally-installed Pixi R, despite it being on the PATH. 
It will not recognize Pixi R even if you point it directly to it. 

Conda environment detections is [unofficially supported](https://github.com/posit-dev/positron/pull/6988), and it works! Put the following in your Positron settings:
```{json}
"positron.r.interpreters.condaDiscovery": true
```

# Known issues trying to use system-wide `R` on `tononi-2`. 
 - `languageserver` should already be installed system-wide, but only Dan can update it!
    - If you install to a personal library, you need to make sure that `.libPaths()` gives priority to the personal library! This can cause other issues. 
    - VSCode will basically never find `languageserver`, even if it's installed right into the system-wide site library and that's first on `.libPaths()`. 
- Installing `httpgd` is a pain. It is usually not available on CRAN. Installing from GitHub in an environment where you don't have root often fails. 
    Sometimes you can go to CRAN and just download the latest version that hasn't been flagged and removed. 
    ```{r}
    install.packages("remotes")
    remotes::install_version("httpgd", "2.0.4")
    ```