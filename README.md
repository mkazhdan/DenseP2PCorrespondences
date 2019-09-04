<center><h1>Dense Point-to-Point Correspondences for Genus-Zero Surfaces (V1.0)</h1></center>
<center>
<a href="#LINKS">links</a>
<a href="#EXECUTABLE">executables</a>
<a href="#USAGE">usage</a>
<a href="#NOTES">notes</a>
<a href="#CHANGES">changes</a>
</center>
<hr>
This distribution contains code for constructing dense point-to-point correspondences between two genus-zero models. Specifically, it provides implementations for:
<ul>
<li>Computing a M&ouml;bius-centered conformal parametrization over the sphere</li>
<li>Evolving the parametrization to make it authalic</li>
<li>Computing Heat Kernel Signatures for a surface</LI>
<li>Rotationally aligning two spherical parametrizations and then refining using optical flow</li>
<li>Transforming a pair of spherical parametrizations into a pair of vertex-to-barycentric-coordinates correspondences</li>
</ul>
<hr>
<a name="LINKS"><b>LINKS</b></a><br>
<ul>
<b>Papers</b>: <a href="http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf">[Kazhdan, Solomon, and Ben-Chen, 2012]</a> <a href="http://www.cs.jhu.edu/~misha/MyPapers/SGP18.pdf">[Baden, Crane, and Kazhdan, 2018]</a> <a href="http://www.cs.jhu.edu/~misha/MyPapers/SGP19.pdf">[Lee and Kazhdan, 2019]</a><br>
<b>Binaries</b>: <a href="http://www.cs.jhu.edu/~misha/Code/DenseP2PCorrespondences/Version1.0/DenseP2PCorrespondences.x64.zip">Windows Executables</a> <a href="http://www.cs.jhu.edu/~misha/Code/DenseP2PCorrespondences/Version1.0/DenseP2PCorrespondences.x64.lib.zip">Windows libraries and DLLs</a><br>
<b>Source Code</b>: <a href="http://www.cs.jhu.edu/~misha/Code/DenseP2PCorrespondences/Version1.0/DenseP2PCorrespondences.zip">ZIP</a> <a href="https://github.com/mkazhdan/DenseP2PCorrespondences">GitHub</a><br>
<b>License</b>: <a href="http://www.cs.jhu.edu/~misha/Code/DenseP2PCorrespondences/license.txt">BSD</a><br>
</ul>

<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>

<ul>
<dl>
<details>
<summary>
<font size="+1"><b>CMCFViewer</b></font>:
Computes a conformal parametrization of a water-tight, genus-zero surface to the sphere using Conformalized Mean Curvature Flow <a href="http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf">[Kazhdan, Solomon, and Ben-Chen, 2012]</a> and  cannonically centers the parametrization relative to M&ouml;bius inversions <a href="http://www.cs.jhu.edu/~misha/MyPapers/SGP18.pdf">[Baden, Crane, and Kazhdan, 2018]</a>. If no output is specified and OpenGL is supported (i.e. the "NO_OPEN_GL" is not defined during the compilation) the executable opens a viewer showing the evolution of the Conformalized Mean Curvature Flow.
</summary>
<dt><b>--in</b> &lt;<i>input mesh</i>&gt;</dt>
<dd> This string is the name of the file from which the mesh will be read.<br>
The file is assumed to be in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions), "px", "py", "pz" (for the associated positions on the unit sphere), and "red", "green", "blue" (for the per-vertex colors). If no colors are provided, they will be assigned using the surface normals.
</dd>

<dt>[<b>--out</b> &lt;<i>output spherical parametrization</i>&gt;]</dt>
<dd> This string is the name of the file to which the spherically parametrized mesh will be written.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions), "px", "py", "pz" (for the associated positions on the unit sphere), and "red", "green", "blue" (for the per-vertex colors).
</dd>

<dt>[<b>--steps</b> &lt;<i>number of CMCF iterations</i>&gt;]</dt>
<dd> This integer values specifies the number of Conformalized Mean Curvature Flow iterations to be used to obtain the conformal spherical parametrization.<br>
The default value for this parameter is 100.
</dd>

<dt>[<b>--stepSize</b> &lt;<i>the temporal size of each CMCF step</i>&gt;]</dt>
<dd> This floating point values specifies the units for the temporal discretization of the Conformalized Mean Curvature Flow.<br>
The default value for this parameter is 0.1.
</dd>

<dt>[<b>--cutOff</b> &lt;<i>Möbius centering cut-off</i>&gt;]</dt>
<dd> This floating point value specifies the threshold for terminating the Möbius centering iterations.<br>
The default value for this parameter is 10^(-10).
</dd>

<dt>[<b>--c2i</b> &lt;<i>center to inversion type</i>&gt;]</dt>
<dd> This integer value specifies how the gradient descent value is to be interpreted as a center of inversion. A value of <b>0</b> indicates that a trivial interpretation is to be used. A value of <b>1</b> indicates that a golden section search should be performed along the descent direction. A value of <b>2</b> indicates that the length of centering transformation should be rescaled using the metric for the Poincar&eacute; disk model.<br>
The default value for this parameter is 2.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If enabled, details regarding the running times of the different stages of processing are output.
</dd>

</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>AuthalicEvolutionViewer</b></font>:
Authalically evolves the parametrization so that mass is uniformly distributed over the sphere and performs a few steps of M&ouml;bius centering at each iteration. If no output is specified and OpenGL is supported (i.e. the "NO_OPEN_GL" is not defined during the compilation) the executable opens a viewer showing the authalic evolution.
</summary>
<dt><b>--in</b> &lt;<i>input spherical parametrization</i>&gt;</dt>
<dd> This string is the name of the file from which the spherically parametrized mesh will be read.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions), "px", "py", "pz" (for the associated positions on the unit sphere), and "red", "green", "blue" (for the per-vertex colors). If the input contains colors, they will be copied to the output. Otherwise, colors are assigned using the surface normals.
</dd>

<dt>[<b>--out</b> &lt;<i>output spherical parametrization</i>&gt;]</dt>
<dd> This string is the name of the file to which the authalically evolved spherical parametrization mesh will be written.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions), "px", "py", "pz" (for the associated positions on the unit sphere), and "red", "green", "blue" (for the per-vertex colors). If the input contains colors, they will be copied to the output. Otherwise, colors are assigned using the surface normals.
</dd>

<dt>[<b>--res</b> &lt;<i>spherical grid resolution</i>&gt;]</dt>
<dd> This integer values specifies the resolution of the spherical grid over which the parametrization is discretized.<br>
The default value for this parameter is 256.
</dd>

<dt>[<b>--steps</b> &lt;<i>number of authalic evolution iterations</i>&gt;]</dt>
<dd> This integer values specifies the number of evolution iterations to be used to make the distribution of mass over the sphere become uniform.<br>
The default value for this parameter is 100.
</dd>

<dt>[<b>--stepSize</b> &lt;<i>temporal size of each authalic evolution step</i>&gt;]</dt>
<dd> This floating point values specifies the units for the temporal discretization of the authalic evolution.<br>
The default value for this parameter is 0.005.
</dd>

<dt>[<b>--cIters</b> &lt;<i>number of centering iterations to be performed in each step</i>&gt;]</dt>
<dd> This integer values specifies the number of M&ouml;bius centering steps that are to be performed after each step of authalic evolution.<br>
The default value for this parameter is 3.
</dd>

<dt>[<b>--smooth</b> &lt;<i>heat diffusion value</i>&gt;]</dt>
<dd> This floating point values specifies the amount of smoothing to applied to the log scale factors prior to computing gradients.<br>
The default value for this parameter is 0.0025.
</dd>

<dt>[<b>--useGSS</b>]</dt>
<dd> If enabled, golden section search should be performed for M&ouml;bius centering. Otherwise, the Poincar&eacute; disk model is used.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If enabled, details regarding the running times of the different stages of processing are output.
</dd>

</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>GetHKS</b></font>:
Performs fast spherical correlation to align the rotational component of the Möbius transformation and outputs the transformed source. (The transformed source is obtained by rotating the parametric positions and, optionally, by replacing the source vertex positions with the corresponding target vertex positions.)
</summary>
<dt><b>--in</b> &lt;<i>input mesh</i>&gt;</dt>
<dd> This string is the name of the file from which the mesh will be read.<br>
The file is assumed to be in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format.<br>
</dd>

<dt>[<b>--out</b> &lt;<i>output Heat Kernel Signatures</i>&gt;]</dt>
<dd> This string is the name of the file to which the Heat Kernel Signatures will be written.<br>
</dd>

<dt>[<b>--dim</b> &lt;<i>number of Heat Kernel Signatures time scales</i>&gt; &lt;<i>the values of the HKS time scales</i>&gt;]</dt>
<dd> This collection of integer and floating points specifies the number of Heat Kernel Signatures to be computed and the time scale for each one.<br>
The default value for this parameter size is 6 Heat Kernel Signatures with time scales 1., 0.5, 0.25, 0.125, 0.0725, and 0.01.<br>
</dd>

<dt>[<b>--dim</b> &lt;<i>spectrum dimension</i>&gt;]</dt>
<dd> This integer value specifies the maximum number of eigenvectors to be used in computing the Heat Kernel Signature.<br>
The default value for this parameter is 200.<br>
</dd>

<dt>[<b>--off</b> &lt;<i>spectral offset</i>&gt;]</dt>
<dd> This floating point value specifies the shift used in the shift-and-inverse implementation of the generalized eigenvalue problem.<br>
The default value for this parameter is 100.0.<br>
</dd>

<dt>[<b>--lump</b>]</dt>
<dd> If enabled, the mass matrix used for computing the spectrum is lumped..<br>
[This is only used in the case that the input is in PLY format.]
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If enabled, details regarding the running times of the different stages of processing are output.
</dd>

</details>
</dl>
</ul>

<ul>
<dl>
<details>
<summary>
<font size="+1"><b>OpticalFlowViewer</b></font>:
Performs the rotational alignment and optical flow refinement to align two spherical parametrizations. If no output is specified and OpenGL is supported (i.e. the "NO_OPEN_GL" is not defined during the compilation) the executable opens a viewer showing the optical flow.<br>
In addition to source and target spherical parametrizations, the executable requires source and target signals. Assigning <I>L</I>+2 values to each vertex. The first value represents the signal used for rotational alignment. The next <I>L</I> are used for optical flow, with one signal being used at each of the <I>L</I> levels of the multiresolution hierarchy. The final one is used for selecting between the rotational candidates.
</summary>
<dt><b>--in</b> &lt;<i>source and target spherical parametrizations</i>&gt;</dt>
<dd> This pair of strings give the names of the files from which the source and target spherically parametrized meshes will be read.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions), "px", "py", "pz" (for the associated positions on the unit sphere), and "red", "green", "blue" (for the per-vertex colors). If the input contains colors, they will be copied to the output. Otherwise, colors are assigned using the surface normals.
</dd>

<dt>[<b>--signal</b> &lt;<i>source and target signals for registration</i>&gt;]</dt>
<dd> This pair of strings give the names of the files containing the (HKS) signals used for registration.
</dd>

<dt>[<b>--out</b> &lt;<i>output source and target parametrizations</i>&gt;]</dt>
<dd> This pair of strings give the names of the files to which the registered spherical parametrizationd mesh will be written.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions), "px", "py", "pz" (for the associated positions on the unit sphere), and "red", "green", "blue" (for the per-vertex colors).
</dd>

<dt>[<b>--res</b> &lt;<i>initial spherical grid resolution</i>&gt;]</dt>
<dd> This integer values specifies the resolution of the coarsest spherical grid used for registration. At each subsequent level of the multi-resolution hierarchy, the resolution is multiplied by a factor of two.<br>
The default value for this parameter is 16.
</dd>

<dt>[<b>--rotRes</b> &lt;<i>spherical grid resolution for rotational alignment</i>&gt;]</dt>
<dd> This integer values specifies the resolution of the spherical grid used for rotational alignment. At each subsequent level of the multi-resolution hierarchy, the resolution is multiplied by a factor of two.<br>
If this parameter is not specified, the finest resolution of the multiresolution hierarchy is used.
</dd>

<dt>[<b>--sWeight</b> &lt;<i>flow field smoothing weight</i>&gt;]</dt>
<dd> This floating point value value specifies given to the flow field smoothness term in the optical flow.<br>
The default value for this parameter is 0.05.
</dd>

<dt>[<b>--solvesPerLevel</b> &lt;<i>number of optical flow iterations</i>&gt;]</dt>
<dd> This integer values specifies the number of optical flow iterations to be performed at each level of the multiresolution hierarchy..<br>
The default value for this parameter is 6.
</dd>

<dt>[<b>--maxRot</b> &lt;<i>maximum number of rotation candidates</i>&gt;]</dt>
<dd> This integer value specifies the maximum of rotation candidates to consider for alignment.<br>
The default value for this parameter is 4.
</dd>

<dt>[<b>--rSeparation</b> &lt;<i>rotation separation distance</i>&gt;]</dt>
<dd> This floating point value specifies the minimum Frobenius distance between rotations when computing rotational alignment candidates.<br>
The default value for this parameter is 3.0.
</dd>

<dt>[<b>--cFraction</b> &lt;<i>minimum correlation fraction</i>&gt;]</dt>
<dd> This floating point value specifies the lower bound on fraction between the correlation value of a candidate rotation and the maximal corelation value.<br>
The default value for this parameter is 0.9.
</dd>

<dt>[<b>--rotate</b> &lt;<i>rotation candidate index</i>&gt;]</dt>
<dd> This integer value specifies which rotation candidate to use for registration.<br>
If no value is specified and output files are specified, all rotation candidates are considered and the result from the best one are written out.
If no value is specified and the viewer is used, no rotation is applied.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If enabled, details regarding the running times of the different stages of processing are output.
</dd>

</details>
</dl>
</ul>

<ul>
<dl>
<details>
<summary>
<font size="+1"><b>GetCorrespondences</b></font>:
Transforms registered spherical parametrizations to vertex correspondences.
</summary>
<dt><b>--in</b> &lt;<i>registered source and target spherical parametrizations</i>&gt;</dt>
<dd> This pair of strings give the names of the files from which the registered source and target spherically parametrized meshes will be read.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions) and "px", "py", "pz" (for the associated positions on the unit sphere).
</dd>

<dt>[<b>--out</b> &lt;<i>output correspondences</i>&gt;]</dt>
<dd> This pair of strings give the names of the files to which the source-to-target and target-to-source correspondences are written. The source-to-target (resp. target-to-source) correspondence files is written in ASCII with a single line for each source (resp. target) vertex, giving the integer index of the target (resp. source) triangle containing the corresponding point and the three floating point barycentric coordinates.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If enabled, details regarding the running times of the different stages of processing are output.
</dd>

</details>
</dl>
</ul>

<hr>
<a name="USAGE"><b>USAGE</b></a><br>
Given genus-zero source and target mesh files, <i>source.ply</i> and <i>target.ply</i>, correspondences can be computed by running the following steps:
<OL>
<LI>Computing the centered spherical parametrizations:
<blockquote><code>
% CMCFViewer --in source.ply --out source.cmcf.ply
</code></blockquote>
<blockquote><code>
% CMCFViewer --in target.ply --out target.cmcf.ply
</code></blockquote>
<LI>Performing the authalic evolution:
<blockquote><code>
% AuthalicEvolutionViewer --in source.cmcf.ply --out source.cmcf.ae.ply
</code></blockquote>
<blockquote><code>
% AuthalicEvolutionViewer --in target.cmcf.ply --out target.cmcf.ae.ply
</code></blockquote>
<LI>Computing the Heat Kernel Signatures to use as registration signals:
<blockquote><code>
% GetHKS --in source.ply --out source.signal
</code></blockquote>
<blockquote><code>
% GetHKS --in target.ply --out target.signal
</code></blockquote>
<LI>Registering the source and target:
<blockquote><code>
% OpticalFlowViewer --in source.cmcf.ae.ply target.cmcf.ae.ply --signal source.signal target.signal --out source_to_target.ply target_to_source.ply
</code></blockquote>
<LI>Computing the correspondences from the registered spherical parametrizations:
<blockquote><code>
% GetCorrespondences --in source_to_target.ply target_to_source.ply --out source_to_target.txt target_to_source.txt
</code></blockquote>
</OL>

<hr>
<a name="NOTES"><b>NOTES</b></a><br>
<ul>
<li> The implementation of this code relies on the <a href="http://eigen.tuxfamily.org/">Eigen</a>, <a href="https://spectralib.org/">Spectra</A>, <a href="https://www.cs.dartmouth.edu/~geelong/soft/">SOFT</a>, and <a href="http://www.fftw.org/">FFTW</a> libraries. The source for Eigen, Spectra, and SOFT are included and should compile under both Windows and Linux. For the FFTW, Windows .lib and .dll files can be found <a href="http://www.cs.jhu.edu/~misha/Code/DenseP2PCorrespondences/Version1.0/DenseP2PCorrespondences.x64.lib.zip">here</a>.</li></ul>

<hr>
<details>
<summary>
<a name="CHANGES"><b>CHANGES</b></a><br>
</summary>
</details>
