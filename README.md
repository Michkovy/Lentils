# How to use LENTILS

LENTILS is run by executing main.py. After execution Python will ask the user for input data, there are three ways
of entering data:

-   text file
-   manual input
-   batch mode

Input via text file is the most common way to enter data. This can be
done by entering a path to a correctly formated ```.txt``` file. See below for an example. If the file is 
valid, LENTILS will proceed and compute the model parameter.

If no path is entered, the file doesn't exist, or is not
formatted properly the program switches to manual input mode. In this
mode the user can enter the same data as found in the text file on a line by
line basis into the console.

Finally, batch mode is activated if the path
entered leads to a folder instead of a file, in this case LENTILS attempts to proces each ```.txt``` file in the specified folder.

In all cases the data needs to be formatted as shown below. Manual input
mode only differs from file mode in that each line has to be entered
individually, as the syntax is completely the same in both cases.

## Text file input

The required lines in the input file are: Name, Lens coordinates, Image
A properties, Image B properties, etc. Lines starting with \# are
considered comments. White spaces between the input
values are ignored.

Example input files can be found in ```./data/```. Detailed descriptions for each line can be found in the manual input
section.

## Batch Processing

In case the entered path leads to a directory LENTILS will produce a
list of text files within the specified directory and proceed to call
its ```main()``` function with each file as input. This is essentially an
automated single file input mode with the printout suppressed. Results
are saved in the programs root folder in the csv formatted ```results.txt```
See below for information on the formatinglayout. Note that this feature is still in early alpha.
There are no provisions for LENTILS encountering a wrongly formatted
file and how to deal with such an exception. At the moment LENTILS will
simply exit batch mode if such a file is encountered. Further
development plans include increasing the robustness of the batch mode
and a way to read single file multi lens inputs.

## Manual Input

If the data input query results in a string LENTILS can't read, either
due to wrong syntax, non-existent file, or is empty the program will
proceed to manual entry mode. The expected single line inputs follow the
file syntax on a line by line basis, pressing ```ENTER``` advances to the
next line.
The first item entered is the system name. It is a string and will be
used to refer to the lens in the output, this may be left empty.

```
SDSSJ1206+4332
```

Next are the lens coordinates. The values from left to right are

-   right ascension, [as],  float
-   right ascension error, [as], float
-   declination, [as], float
-   declination error, [as], float

The important relation for calculating the model are the relative
image-lens positions, as long as this relation is preserved it doesn't
matter if RA/Dec are absolute or relative measurements. Since LENTILS
centres the system on the given lens position for the computation
internally any of the objects can be chosen as the origin of the
coordinate system on entry.

```
-0.664,.0137,1.748,.028
```

The third line is the start of the image data input. Each image's
coordinates together with its magnitude are entered line by line. It
follows the same "value, value error" pair syntax as the lens data. The
required image inputs are:

-   designation, string
-   right ascension, [as], float
-   right ascension error, [as], float
-   declination, [as], float
-   declination error, [as], float
-   magnitude, [mag], float
-   magnitude error, [mag], float

Designation will be used to link kappa, gamma, mu values to a
specific image in the output. It is also used in the optional plot, it is recommended
to keep it short and distinct to avoid clutter in the plot.
Magnitudes are a required line item for both quads and doubles, although
for quads it can be safely set to zero in v1.0 since they are only used
for determining parameters in double lenses.

Since the magnitude difference is the crucial parameter in the computation, either
apparent magnitude or magnitude differences ($\Delta m$) between the
images work as input. LENTILS assumes $\Delta m=m_B-m_A$ internally, which means when using only magnitude difference, the right way to enter is $m_A=0$
and $m_B=\Delta m$, if $\Delta m$ was computed via $\Delta m=m_B-m_A$.
Flux ratios must be manually converted to magnitude differences using
$\Delta m=-2.5log(f_b/f_a)$.

```
A,0,0.11,0,0.01,18.05,.02
```

LENTILS requires at least 2 image lines. In the case of a double lens
press enter twice, effectively entering empty image coordinates, after
the second image line to start computation, for quads simply enter two
more image lines after the second image line and the program will start
computation after the fourth line is entered.

```
B,-0.098,0.006,2.894,.009,18.38,.02
```

Lens and image data are then passed to the appropriate SIS+$\gamma$
solver, ```Twin2()``` or ```Quad2()``` respectively, which will return:

-   Einstein radius/characteristic scale of the lens R_0
-   Source position relative to lens position vector beta 
-   External shear vector gamma_{ext}

Per default the parameters are calculated thrice, first to find the
parameter itself followed by a second and third run using the
input plus/minus uncertainty as the input for the parameter error calculation
Per convention we assume the upper uncertainty increases
the parameter value. 
If **all** input uncertainties are zero,
LENTILS runs only the parameter computation and sets the ```no_error$```
flag to ```True``` whilst returning the uncertainties as zero. Parameters
along with their uncertainties are saved internally in the ```sis```
dictionary.








The ```lensplot()``` produces of the macromodel displays, input
image position as black $\times$, the lens position as a white cross,
the model derived source position is shown as a blue $\times$, model
predicted image positions are shown similarly as a red cross. Image
positions are the solutions to equation
[\[eq41\]](#eq41){reference-type="ref" reference="eq41"} using the
computed model parameters. Since double lenses produce an exact result
the model images are at the same positions as the observed images, while
quad images show a difference from their observed counterparts. For
every position with an uncertainty there is also an error ellipse
centred on its location. The area around the lens is shown in a square
$2\times125\%$ the size of the largest lens-image separation. Over that
region LENTILS calculates $\mu$
(eq.[\[eq27\]](#eq27){reference-type="ref" reference="eq27"}) on
250$\times$`<!-- -->`{=html}250 grid centred on the lens position. These
two parameters may be changed if the plot routine is called manually,
see A.[\[A_lensplot\]](#A_lensplot){reference-type="ref"
reference="A_lensplot"}, but in practice we've found these values to
provide a good compromise between speed and smooth contours. The
background then shows the filled contours of $log(|(1/\mu)|)$ where
darker areas indicate less magnification. Of particular interest are
locations where $1/\mu$ is infinite since they imply det$(\hat{A})=0$,
which by definition is the critical line demarcating the border between
strong and weak lensing regions. The isopleth det($\hat{A})=0$ is drawn
as a solid red line. We also plot a stroked blue circle with radius
$R_0$ centred on the lens position indicating the Einstein radius of the
lens on the source plane[^1].

While $\kappa$ and $\gamma$ by default are only computed at the image
position, the SISg (A.[\[A_SISgFunc\]](#A_SISgFunc){reference-type="ref"
reference="A_SISgFunc"}) routine is capable of taking in an array of
coordinates so that the parameters may be calculated over the whole
region. The returned variables are scalar fields for convergence, or
total shear magnitude, that can be plotted the same way as $\mu$.

[^1]: By convention image plane objects are red with solid stroke,
    source plane objects are blue and dashed, inputs are black or white
