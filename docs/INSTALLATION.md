Installation {#installpage}
============

### Dependencies

* [htslib](https://github.com/samtools/htslib)
* [Awk](https://www.cs.princeton.edu/~bwk/btl.mirror/) - Available on most UNIX systems.

Note:
* It is highly recommended that [samtools](https://github.com/samtools/samtools) also be installed alongside iVar. iVar uses the output of samtools mpileup to call variants and generate consensus sequences. In addition, samtools `sort` and `index` commands are very useful to setup a pipeline using iVar.


### Installing on Mac

To install ivar, run the following commands.

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

Also, please add the following to your .bash_profile so that iVar can find htslib dynamic libraries during runtime.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/hts/lib/folder
```

Installing [GNU Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html#Autotools-Introduction),


Using [homebrew](https://brew.sh/),

```
brew install autoconf automake libtool
```

### Installing on Linux

To install ivar, run the following commands.

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

Also, please add the following to your .bash_rc so that iVar can find htslib dynamic libraries during runtime.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/hts/lib/folder
```

Installing [GNU Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html#Autotools-Introduction),

Using [APT](https://help.ubuntu.com/lts/serverguide/apt.html) on Ubuntu,

```
apt-get install autotools-dev
```

For bug reports please email gkarthik[at]scripps.edu or raise an issue on Github.
