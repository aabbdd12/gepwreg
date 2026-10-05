# gepwreg 1.5.0 -- the copy submitted to the SSC archive

The same files as gepwreg 1.5.0 on GitHub (tag v1.5.0), flat, with the paper
as an ancillary file (gepwreg_paper.pdf, version 5 of the Zenodo record
10.5281/zenodo.22845707). It is here to test the SSC package before the
archive publishes it. As on SSC, the package file lists every file with an
`f` line: `net install` installs the programs, the help and the dialog only;
the example data and the paper are copied by `net get` to the current folder,
never to the system directories.

```stata
net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/main/ssc/1.5.0") replace
net get     gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/main/ssc/1.5.0") replace
```

The examples of the help read the data from the current folder, else from the
SSC archive, else from GitHub. To go back to the GitHub version:

```stata
net install gepwreg, from("https://raw.githubusercontent.com/aabbdd12/gepwreg/main") replace
```
