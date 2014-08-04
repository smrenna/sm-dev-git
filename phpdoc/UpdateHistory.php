<html>
<head>
<title>Update History</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='UpdateHistory.php'>
 
<h2>Update History</h2> 
 
These update notes describe major updates relative to the 
PYTHIA 8.186 version, which was the last regular 8.1 release.
(Minor bug fixes will continue to appear.) The step from 
8.1 to 8.2 gave an occasion to break backwards compatbility,
but this should only affect a small part of the user code.

<h3>Main news by version</h3> 
 
<ul> 
 
<li>8.200: 4 August 2014 
<ul> 
 
<li>A new <code>pdfdoc</code> directory collects pdf documents
that are linked from the <code>htmldoc</code> and <code>phpdoc</code>
directories. Over time it will  provide more in-depth descriptions of 
various physics aspects than offered in the html/php-formatted 
documentation.</li> 
 
<li>A new <code>include/Pythia8Plugins</code> directory collects
code that does not form part of the core PYTHIA functionality but
still has a general usefulness. Code in this directory will not be
compiled as part of the Pythia library, but can be linked where needed.
This new directory currently contains
<ul>
<li>the jet matching classes in <code>CombineMatchingInput.h</code>,
<code>GeneratorInput.h</code> and <code>JetMatching.h</code>, moved
from the <code>examples</code> directory;</li> 
<li>the <code>PowhegHooks</code> user hook, to veto shower emissions 
above the POWHEG scale, formerly found in <code>examples/main31.cc</code>;
</li>
<li>the <code>Pythia8ToHepMC</code> interfoace for output of PYTHIA events 
into the HepMC format, combining the code previously in 
<code>include/Pythia8ToHepMC.h</code> and
<code>pythia8tohepmc/Pythia8ToHepMC.cc</code>; and</li>
<li>the <code>FastJet3.h</code> interface of PYTHIA particles to the
FastJet 3 library of jet finders, formerly found in 
<code>include/FastJet3.h</code>.</li> 
</ul>
 
<li>Several methods have been removed from the <code>Event</code> class
since the properties now instead can be accessed from the individual 
<code>Particle</code> instance, if this particle belongs to an event.
These include <code>iTopCopy</code>, <code>iBotCopy</code>,
<code>iTopCopyId</code>, <code>iBotCopyId</code>,<code>motherList</code>,  
<code>daughterList</code>, <code>sisterList</code>,
<code>sisterListTopBot</code>, <code>isAncestor</code>,
<code>statusHepMC</code> and <code>undoDecay</code>.</li> 
 
<li>A number of deprecated <code>Pythia::init(...)</code> methods with
varying arguments have been removed. Instead call <code>init()</code>
without any arguments and use 
<?php $filepath = $_GET["filepath"];
echo "<a href='BeamParameters.php?filepath=".$filepath."' target='page'>";?>Beam Parameters</a> settings to
specify beams and energies in different ways.</li> 
 
<li>The deprecated <code>Pythia::statistics(...)</code> method has been
removed; instead use <code>Pythia::stat(...)</code>.</li> 
 
<li>Several settings in the <code>Main:</code> series have been removed.
Most of these have already found replacements in the <code>Init:</code>,
<code>Next:</code> and <code>Stat:</code> ones, and have been marked as 
deprecated. Four further ones were deemed so peripheral that they were
removed altogether, but of course the underlying functionality remains.
</li> 
 
<li>A few aliases for (parts of) settings names have been removed.
Previously "Multiple" was mapped to "Multiparton", "MI" to "MPI" and
"minBias" to "nonDiffractive" if a settings name was not found for the 
original inpout string.</li> 
 
<li>The default tune has been changed from 4C to Monash 2013, meaning
<code>Tune:ee = 7</code> and <code>Tune:pp = 14</code>. The old 4C
tune that was default in 8.1 can be recovered with 
<code>Tune:ee = 3</code> and <code>Tune:pp = 5</code>.
Also most other older tunes are based on <code>Tune:ee = 3</code>.
</li> 
 
<li>Two new CMS underlying-event tunes [<a href="Bibliography.php" target="page">CMS14</a>] have been added 
as options.</li> 
 
<li>Christine Rasmussen joins as new PYTHIA collaboration member.</li> 
 
<li>A new model for the handling of <?php $filepath = $_GET["filepath"];
echo "<a href='BeamRemnants.php?filepath=".$filepath."' target='page'>";?>beam 
remnants</a> as an option to the old one, which remains as default
for now.</li> 
 
<li>Two new models for <?php $filepath = $_GET["filepath"];
echo "<a href='ColourReconnection.php?filepath=".$filepath."' target='page'>";?>colour 
reconnection</a>, one quite sophisticated and one simpler.
This involves several new classes and files. It also includes some 
changes in the hadronization framework, notably for the handling of 
junctions. The old model remains as default for now. The 
<code>BeamRemnants:reconnectColours</code> flag to switch on/off
reconnection has been renamed <code>ColourReconnection:reconnect</code>,
the main parameter <code>BeamRemnants:reconnectRange</code> of the old
model has been renamed <code>ColourReconnection:range</code>, and several
new settigns have been introduced, notably 
<code>ColourReconnection:mode</code> to switch among the three models.
</li> 

<li>A new <code>include/Pythia8Plugins/ColourReconnectionHooks.h</code> 
makes available an even larger selection of toy colour reconnection 
models, via user hooks. Some of them are only intended for top decays, 
for top mass uncertainty studies, whereas others can be used more 
generally. The <code>examples/main29.cc</code> program illustrated how 
the different options should be set up.</li> 

<li>Improved capability for the <code>LHAup</code> Les Houches interface 
to read SLHA information embedded in the input file or stream.</li> 

<li>The <code>Makefile</code>s have been updated to take into account
the changed structure of the HepMC interface.</li> 
 
<li>Bug fix in the two-loop running <i>alpha_s</i>, for the matching 
to six flavours at the top mass.</li> 
 
<li></li> 
 
</ul> 
</li> 
 

</ul> 
 
</body>
</html>
 
<!-- Copyright (C) 2014 Torbjorn Sjostrand --> 
