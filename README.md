iVar
===========

iVar is a computational package that contains functions broadly useful for viral amplicon-based sequencing. Additional tools for metagenomic sequencing are actively being incorporated into iVar. While each of these functions can be accomplished using existing tools, iVar contains an intersection of functionality from multiple tools that are required to call iSNVs and consensus sequences from viral sequencing data across multiple replicates. We implemented the following functions in iVar: (1) trimming of primers and low-quality bases, (2) consensus calling, (3) variant calling - both iSNVs and insertions/deletions, and (4) identifying mismatches to primer sequences and excluding the corresponding reads from alignment files.

## [Manual](https://andersen-lab.github.io/ivar/html/)

## Installation

### Dependencies

* [htslib](http://www.htslib.org/download/)
* [Awk](https://www.cs.princeton.edu/~bwk/btl.mirror/) - Pre-installed on most UNIX systems.

Note:
* It is highly recommended that [samtools](https://github.com/samtools/samtools) also be installed alongside iVar. iVar uses the output of samtools mpileup to call variants and generate consensus sequences. In addition, samtools `sort` and `index` commands are very useful to setup a pipeline using iVar.


### Installing on Mac

[GNU Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html#Autotools-Introduction) is required to compile iVar from source.

To install Autotools using [homebrew](https://brew.sh/) please use the command below,

```
brew install autoconf automake libtool
```

To install iVar, run the following commands.

```
./autogen.sh
./configure
make
make install
```

If htslib has been installed in a non standard location, please run the following commands,

```
./autogen.sh
./configure --with-hts=/prefix/to/bin/folder/with/htslib
make
make install
```

If htslib is installed in a non standard location, please add the following to your .bash_profile so that iVar can find htslib dynamic libraries during runtime.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/hts/lib/folder
```

To test installation just run, `ivar version` and you should get the following output,

```
iVar version 1.0

Please raise issues and bug reports at https://github.com/andersen-lab/ivar/
```

### Installing on Linux

[GNU Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html#Autotools-Introduction) is required to compile iVar from source.

To install Autotools using [APT](https://help.ubuntu.com/lts/serverguide/apt.html) please use the command below,

```
apt-get install autotools-dev
```

To install iVar, run the following commands.

```
./autogen.sh
./configure
make
make install
```

If htslib has been installed in a non standard location, please run,

```
./autogen.sh
./configure --with-hts=/prefix/to/bin/folder/with/htslib
make
make install
```

If htslib is installed in a non standard location, please add the following to your .bashrc so that iVar can find htslib dynamic libraries during runtime.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/hts/lib/folder
```

To test installation just run, `ivar version` and you should get the following output,

```
iVar version 1.0

Please raise issues and bug reports at https://github.com/andersen-lab/ivar/
```

### Contact

For bug reports please email gkarthik[at]scripps.edu or raise an issue on Github.
