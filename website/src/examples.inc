<body>
<div id="wrap">
  <div id="top">
    <h2> <a href="index.html"><strong>VCFtools</strong></a></h2>
    <div id="menu">
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><a href="http://sourceforge.net/projects/vcftools/">Sourceforge page</a></li>	
        <li><a href="examples.html" class="current">Examples &amp; Documentation</a></li>
        <li><a href="downloads.html">Downloads</a></li>
      </ul>
    </div>
  </div>
  <div id="content">
    <div id="left">
    	<h2>Getting Started with VCFtools</h2>
		<p>To download VCFtools, please see the <a href="downloads.html">downloads</a> page.
		Once downloaded, open a terminal and change the working directory to the file download location.</p>
		<p>If the program was downloaded as a tarball, it must be uncompressed using the following command (replacing it with current version number):</p>
		<p class="codebox">tar -xvf vcfools.0.X.XX.tar.gz</p>
		<p>Certain scripts within VCFtools require that VCF files are compressed by bgzip and indexed by tabix (both tools are part of the tabix package, available for download <a href="https://sourceforge.net/projects/samtools/files/tabix/">here</a>).
		 Both tools must be in directories that are listed in the PATH environment variable. In addition, to run the VCFtools Perl scripts, the PERL5LIB environment variable must be set to include the Vcf.pm module. This can be achieved as follows:</p>
		<p class="codebox">export PERL5LIB=/path/to/your/vcftools-directory/perl/</p>
		<p>The program must now be compiled. Change to the specified target directory.</p>
		<p class="codebox">cd vcftools/
		<br>make</p>
		<p>This will compile both the PERL API and the C++ executable. By default, the compiled programs will now be in the <b>vcftools/bin/</b> directory.
		Some common compiler errors may appear, such as not being able to locate <b>Vcf.pm</b> or not being able to locate <b>zlib</b>.
		If these occur, please visit the <a href="http://zlib.net/">zlib page</a> for installation instructions or edit the PERL5LIB variable as explained above.</p>

      <h2>Specific examples by module</h2>
      <h3><a href="perl_examples.html">Usage Examples for the PERL API</a></h3>
      <br>
      <h3><a href="documentation.html">Usage Examples for the Binary Executable</a></h3>
		
    </div>
     <div id="right">
      <div class="box">
      	<h2>Binary executable documentation</h2>
      		<p>The manual page for:</p>
      	<p><a href="man_latest.html">v0.1.13</a> - developmental version</p>
		<p><a href="man_0112b.html">v0.1.12b</a> - current released version</p>
		<p><a href="man_0112a.html">v0.1.12a</a></p>
      	<h2>PERL module documentation</h2>
      		<p>The full documentation for the <a href="perl_module.html">latest version</a>.</p>
      	<h2>HTSlib</h2>
      		<p>Documentation for the separate HTSlib package is <a href="htslib.html">here</a>.</p>
      </div>
    </div>
    <div id="clear"></div>
  </div>
  <div id="footer">
    <p>Copyright 2014 &copy; VCFtools</p>
  </div>
</div>
</body>
</html>
