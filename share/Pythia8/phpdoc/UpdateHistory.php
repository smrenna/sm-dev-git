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
8.1 to 8.2 gave an occasion to break backwards compatibility, 
but this should only affect a small part of the user code. 
 
<h3>Main news by version</h3> 
 
<ul> 
 
<li>8.222: 25 December 2016 
<ul> 
 
<li>Nadine Fischer and Leif L&ouml;nnblad join as co-authors, 
while Jesper Roy Christiansen leaves. Nishista Desai, Ilkka Helenius 
and Stefan Prestel have new affiliations.</li> 
 
<li>The machinery for resolved &gamma;&gamma; collisions has been extended, 
such that now soft processes and MPIs can be simulated, also when 
embedded in <i>l^+l^-</i> collisions. (But not yet diffraction.) 
Also some further improvements have been introduced, see 
the <?php $filepath = $_GET["filepath"];
echo "<a href='PhotonPhoton.php?filepath=".$filepath."' target='page'>";?>Photon-photon Interactions</a> 
description. This implies several changes in different parts of 
the code, mainly related to beam remnants and beam particles.</li> 
 
<li> 
Also direct-resolved and direct-direct processes are included for 
&gamma;&gamma; interactions, with photon beams and within lepton beams. 
This involves new subprocesses where one initiator is a photon and the 
other a parton. A new sample main program (<code>main69.cc</code>) 
illustrates how the different classes of &gamma;&gamma; interactions 
are combined.</li> 
 
<li>The kinematics of &gamma;&gamma; have been revised to include all 
mass corrections and to handle also non-equal leptons. A new class 
<code>GammaKinematics</code> is introduced to handle the sampling of 
the kinematics. A fix for the <code>ProcessLevel::roomForRemnants()</code> 
function, which rejected a bit too many processes when photon-photon 
collisions were generated within lepton beams. 
</li> 
 
<li>New cuts for the kinematics of &gamma;&gamma; interactions in 
<i>l^+l^-</i> collisions are introduced, for details see 
<?php $filepath = $_GET["filepath"];
echo "<a href='PhotonPhoton.php?filepath=".$filepath."' target='page'>";?>Photon-photon Interactions</a>. 
Matching new kinematics output methods, see 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Event Information</a>.</li> 
 
<li> 
A PDF for point-like photon is included. In case of lepton PDFs, the 
photon contribution has now restricted virtuality and also more accurate 
lower limit for the virtuality. A new option to use separate PDFs for 
hard processes, with photon PDFs obtained from LHAPDF5. CJKL PDFs are 
modified so that, instead of freezing the scale below its minimum, 
the scale evolution is approximated with <i>log(Q^2)</i>.</li> 
 
<li>A new alternative "thermal hadronization" option is introduced, 
wherein an exponential <i>exp(-pT / T)</i> hadronic transverse 
momentum spectrum replaces the default Gaussian one, with a 
"temperature" <i>T</i> as free parameter. Given this <i>pT</i>, 
the next hadron (consistent with local flavour conservation) is picked 
among the possibilities with an <i>exp(-mT / T)</i> weight. 
This option is accessed with <code>StringPT:thermalModel = on</code>. 
See further the <?php $filepath = $_GET["filepath"];
echo "<a href='Fragmentation.php?filepath=".$filepath."' target='page'>";?>Fragmentation</a> and 
<?php $filepath = $_GET["filepath"];
echo "<a href='FlavourSelection.php?filepath=".$filepath."' target='page'>";?>Flavour Selection</a> descriptions 
and the article [<a href="Bibliography.php" target="page">Fis16</a>].</li> 
 
<li>A new option <code>StringPT:closePacking = on</code> allows to 
enhance the <i>pT</i> width in regions where there is a high 
density of partly overlapping strings. This works both for the 
default Gaussian and the alternative exponential (see above) 
<i>pT</i> description; in the latter case it will also enhance the 
rate of heavier-particle production. See further the 
<?php $filepath = $_GET["filepath"];
echo "<a href='Fragmentation.php?filepath=".$filepath."' target='page'>";?>Fragmentation</a> description 
and the article [<a href="Bibliography.php" target="page">Fis16</a>].</li> 
 
<li>A new simple model for hadronic rescattering is introduced, 
with two variants, as described in [<a href="Bibliography.php" target="page">Fis16</a>]. A new master 
switch <code>HadronLevel:HadronScatter</code>, by default off, and 
<code>HadronScatter:mode</code> to pick among them and the old one. 
See the <?php $filepath = $_GET["filepath"];
echo "<a href='HadronScattering.php?filepath=".$filepath."' target='page'>";?>Hadron Scattering</a> page 
for further details and parameters. New status codes 111 and 112 
introduced for rescattered hadrons.</li> 
 
<li>A first process for Dark Matter production, for the pair production 
of new Dirac DM particles through an <i>s</i>-channel vector-like 
mediator.</li> 
 
<li>New mode <code>BeamRemnants:unresolvedHadron</code> can be used to 
let original hadron remain as beam remnant, e.g. for coherent emission 
of photons.</li> 
 
<li>New option with running coupling in Hidden Valley scenarios. 
Some other small fixes in it.</li> 
 
<li>New behaviour of <code>Event::remove</code>, where mother and 
daughter indices now are updated by default.</li> 
 
<li>Fix in the setup of tunes with 
<code>SpaceShower:rapidityOrder = off</code>. 
The new (in 8.219) <code>SpaceShower:rapidityOrderMPI</code> then also 
ought to have been set off, but this was missed, giving small 
inconsistencies (around 2% reduction of the charged multiplicity). 
Thanks to James Monk.</li> 
 
<li>New method <code>string Hist::getTitle()</code> to get the title 
of a histogram, while <code>void Hist::title(string )</code> sets it. 
Thanks to Roberto Franceschini.</li> 
 
<li>Corrected behaviour for <i>R</i>-hadrons produced in sequential 
resonance decays (for example a squark decaying to a gluino with the 
latter forming an <i>R</i>-hadron). Thanks to Jinmian Li for 
alerting us. 
 
<li>Minor updates so that <code>main91</code> and <code>main92</code> 
examples work also with ROOT 6, in addition to the existing ROOT 5 
support. Thanks to Li Huang.</li> 
 
<li>Include correct mass suppression for excited fermion three-body 
decays <i>F^* &rarr; F Fbar f</i>, primarily for <i>F = t</i>. 
Thanks to Olya Igonkina and Oleg Zenin.</li> 
 
<li>Two new <code>Hist::rivetTable</code> methods allow histograms to 
be written on file in a format that Rivet understands.</li> 
 
<li>New particle data method <code>nQuarksInCode(int idQ)</code> counts 
how many copies of the requested quark code <code>idQ</code> that a 
quark, diquark, meson or baryon code contains.</li> 
 
<li>Bring the <code>FJcore</code> package inside the <code>Pythia8</code> 
namespace to avoid potential name clashes with user code. Thanks to 
Andy Buckley.</li> 
 
<li>New <code>#define PYTHIA_VERSION_INTEGER 82xx</code> in 
<code>Pythia.h</code> matches already existing 
<code>#define PYTHIA_VERSION 8.2xx</code>, for more convenient 
matching using integers. Thanks to Andrii Verbytskyi.</li> 
 
<li>The handling of the <code>meMode</code> ranges 52 - 60 and 62 - 70 
were incorrect, insofar as checks or not against duplication of existing 
channels go, and have now been set straight. Thanks to Christopher 
West.</li> 
 
<li>A minor fix for <code>BeamParticle::popBack()</code> method to reset 
companion choice also if <code>iComp = 0</code>.</li> 
 
<li>Two minor particle data fixes.</li> 
 
<li>Small fix in the parsing code of <code>LHEF3.h</code>.</li> 
 
<li>Small clarifications in the documentation.</li> 
 
</ul> 
</li> 
 
<li>8.219: 10 May 2016 
<ul> 
 
<li>An interface to the Python programming language has been introduced, 
see <?php $filepath = $_GET["filepath"];
echo "<a href='APythonInterface.php?filepath=".$filepath."' target='page'>";?>A Python Interface</a> for details. 
Various minor changes in the C++ code have been done in order to permit 
the automatic generation of the interface.</li> 
 
<li>Included a new framework for automated parton-shower uncertainty bands. 
Variations of the QCD renormalisation scale for both initial- and 
final-state showers can now be computed by Pythia on the fly, and are 
provided as a list of alternative weights for each event, representing 
the probability that the given event would have occurred under different 
shower assumptions. So-called "nonsingular terms" can also be added to 
the splitting kernels to estimate the possible effect of missing 
matrix-element corrections. Full documentation is found on the 
<?php $filepath = $_GET["filepath"];
echo "<a href='Variations.php?filepath=".$filepath."' target='page'>";?>Automated Shower Variations</a> page, 
and a paper is due to appear on arXiv shortly.</li> 
 
<li>When a final-state <i>g &rarr; g g</i> branching happens with 
a massive recoiler, radiation in the recoiler direction is now 
by default further suppressed to respect the "dead cone" effect, 
see new switch <code>TimeShower:recoilDeadCone</code>. Furthermore 
a new switch, <code>TimeShower:MEextended</code>, on by defaults, 
attempts to guess the most relevant ME correction when the 
correct choice is not known or implemented. Thanks to 
Jesse Thaler, Michele Selvaggi and Fabio Maltoni.</li> 
 
<li>The <i>gamma-gamma</i> hard-process machinery has been extended to 
convolute the partonic PDFs in a photon with the flux of photons inside 
a lepton. The description is intended for the region of quasireal photons, 
but full kinematics is implemented.</li> 
 
<li>New constructors that take streams rather than files as input. 
Thus the contents of a file can be read and then broadcast to multiple 
instances of Pythia, eliminating the inefficiency of multiple jobs 
reading from the same initialization files. This facilitates running 
with MPIs (Message Process Interfaces). Affected areas are settings, 
particle data and PDF grids.</li> 
 
<li>Small additions and fixes to the LHEF3 framework, to keep track of 
weight keys, and improve parsing.</li> 
 
<li>A new class allows PDF data files in the lhagrid1 format, 
with some restrictions, to be read and used. This does not replace 
all that you can do with a complete LHAPDF6 installation, but at 
least permits some simple studies without LHAPDF6 + Boost. 
New example <code>main55.cc</code> illustrates this, and event 
properties for an intermediate spinless resonance in 
<i>&gamma; + &gamma; &rarr; &gamma; + &gamma;</i> at 750 GeV.</li> 
 
<li>Four Pomeron PDF sets from the ACTW study [<a href="Bibliography.php" target="page">Alv99</a>] now 
implemented.</li> 
 
<li>The three Pomeron H1 Jets PDF data files have been joined to one. 
</li> 
 
<li>The <?php $filepath = $_GET["filepath"];
echo "<a href='ExternalDecays.php?filepath=".$filepath."' target='page'>";?>external decays</a> interface has 
been extended by a new method that allows sequential decays to be done 
in one call.</li> 
 
<li>Increased possibilities to set the <i>epsilon</i> and 
<i>alpha'</i> parameters of the Pomeron trajectory for 
hard-diffraction Pomeron fluxes.</li> 
 
<li>Extrapolation of PDFs to small <i>x</i> values when 
<code>PDF:extrapolate = on</code> now extended to more cases.</li> 
 
<li>New flag <code>TimeShower:QEDshowerByOther</code> allows charged 
resonances, like the <i>W^+-</i>, to radiate photons.</li> 
 
<li>The jet matching algorithm in <code>JetMatching.h</code> has been 
extended to better handle heavy quarks, heavy colored particles (such 
as squarks) and "other" partons (coloured but produced from an 
Electroweak vertex).</li> 
 
<li>Added a new option <code>SpaceShower:rapidityOrderMPI</code>, which 
will enforce an ordering in rapidity for emissions off secondary 
scattering systems. This will be enabled by default, ensuring backwards 
compatibility. The old <code>SpaceShower:rapidityOrder</code> now only 
refers to the hard(est) subprocess.</li> 
 
<li>Changes to the cross section handling in the presence of user 
vetoes/weights, containing three changes: 
<br/>1. Counter of selected event <code>pythia.info.nSelected()</code> 
is now updated  immediately after the hard process generation. 
<br/>2. More fine-grained input settings to enforce that Pythia 
generates (or reads) exactly a fixed number of hard process events. 
<br/>3. The Pythia "internal cross section" 
<code>pythia.info.sigmaGen()</code> and the event weight 
<code>pythia.info.weight()</code> now directly include the effect of 
event vetoes and event reweighting that are applied by the leading-order 
ME+PS merging prescriptions.</li> 
 
<li>Changes to the merging classes to allow for a postponed CKKW-L 
event veto. (See end of <?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L Merging</a> 
page for details.)</li> 
 
<li>Small changes to virtual parton shower functions, and to the merging 
classes, to prepare ME+PS merging with the Vincia parton shower.</li> 
 
<li>The <?php $filepath = $_GET["filepath"];
echo "<a href='HepMCInterface.php?filepath=".$filepath."' target='page'>";?>HepMC</a> interface has been 
modified such that the detection of unhadronized quarks or gluons 
leads to an exception being thrown, so that the user can decide what 
action to take. See further new/renamed <code>free_parton_exception</code> 
switch. The usage of exceptions for this specific task is by request 
from ATLAS, and does not represent a general change of programming 
style. Thanks to James Monk.</li> 
 
<li>A class <code>WVec</code> has been introduced to store vectors of 
strings. The delimiters { } are introduced to provide input with 
embedded blanks, and broken up across several lines. See further the 
<?php $filepath = $_GET["filepath"];
echo "<a href='SettingsScheme.php?filepath=".$filepath."' target='page'>";?>Settings Scheme</a> description.</li> 
 
<li>The <code>Settings::toLower</code> method used to convert a string 
to lowercase, and also trim it from initial or trailing blanks and special 
characters, now moved to <code>PythiaStdlib.h</code> so it can be used 
more generally. Other code changes accordingly.</li> 
 
<li>Remove many rarely (if ever) used <code>ostream& os = cout</code> 
optional arguments in favour of hardcoded <code>cout</code> in the code. 
Eliminates some redundancy of methods.</li> 
 
<li>Rename <code>...print(...</code> methods to <code>...list(...</code> 
to favour a more regular naming pattern.</li> 
 
<li>Minor <code>configure</code> and <code>Makefile</code> updates, 
to address potential linking problems on some platforms for boost, gzip 
and promc. Thanks to Dmitri Konstantinov.</li> 
 
<li>New <code>--config</code> option of <code>pythia8-config</code> 
echoes the arguments passed to <code>configure</code>.</li> 
 
<li>The description of <code>POWHEG:vetoCount = 0</code> has been 
corrected. Thanks to Florian Koenig.</li> 
 
<li>Allow for lightest neutralino not to decay as a resonance.</li> 
 
<li>Do not switch off the Breit-Wigner width treatment of a resonance 
as easily as previously, but only if the width is below 1e-6 GeV.</li> 
 
<li>Make some <code>SlowJet</code> methods virtual to allow derived 
classes with modified properties.</li> 
 
<li>Moved some misplaced info on parton-level choice of MPIs.</li> 
 
<li>Corrected typos where some bottomonium long-distance matrix element 
had been set larger than normally assumed.</li> 
 
<li>Fixed typo potentially giving incorrect colour flow in resonance 
decays.</li> 
 
<li>Fixed an out-of-bounds array access in HelicityMatrixElements. 
Thanks to Vittorio Zecca.</li> 
 
<li>Fixes in <code>main80.cc</code> and <code>main89.cc</code>.</li> 
 
<li>Fixed problem with the SLHAinterface not being zeroed-out when 
using repeated subruns.</li> 
 
<li>Minor fix for beam particles, that no default was set as to whether 
they are gammas or not.</li> 
 
<li>Cleaned up error printout for the PDF classes.</li> 
 
<li>Change a few true/false to on/off in the documentation to make the 
php version of the manual recognize them.</li> 
 
<li>Some trivial code, manual and bibliography updates.</li> 
 
</ul> 
</li> 
 
<li>8.215: 4 January 2016 
<ul> 
 
<li>Ilkka Helenius joins as new PYTHIA co-author.</li> 
 
<li>A new machinery for <i>gamma-gamma</i> collisions is now available, 
see <?php $filepath = $_GET["filepath"];
echo "<a href='PhotonPhoton.php?filepath=".$filepath."' target='page'>";?>Photon-photon Interactions</a>. 
So far only hard processes can be generated, along with parton showers 
and hadronization, but without multiparton interactions. The CJKL parton 
distributions of the photon have been implemented and are used.</li> 
 
<li>Double production of charmonium and bottomonium <i>3S1</i> states 
is now available, but with only the colour-singlet processes included, 
see <?php $filepath = $_GET["filepath"];
echo "<a href='OniaProcesses.php?filepath=".$filepath."' target='page'>";?>Onia Processes</a> for details.</li> 
 
<li>Weak merging implemented, i.e. <i>W</i> gauge bosons can be 
produced either as part of the hard matrix element or in the parton 
shower, and a proper treatment merges these two possibilities consistently. 
See the <?php $filepath = $_GET["filepath"];
echo "<a href='CKKWLMerging.php?filepath=".$filepath."' target='page'>";?>CKKW-L</a> page for details.</li> 
 
<li>Running <i>alpha_em</i> in merging description.</li> 
 
<li>Improved interface to external parton showers, such as 
<a href="http://vincia.hepforge.org/" target="_top">VINCIA</a> and 
<a href="http://www.slac.stanford.edu/~prestel/DIRE/" target="_top">DIRE</a>, 
so that these now also can use the various matching and merging 
frameworks implemented in Pythia.</li> 
 
<li>New options in the <?php $filepath = $_GET["filepath"];
echo "<a href='JetMatching.php?filepath=".$filepath."' target='page'>";?>jet matching</a> 
framework, such that expert users can use their own veto code for 
Madgraph-style jet matching.</li> 
 
<li>New convenient possibility to run Madgraph5_aMC@NLO from within 
Pythia, by wrapping the Madgraph5_aMC@NLO executable inside a new 
<code>LHAupMadgraph</code> class that derives from the Pythia 
<code>LHAup</code> base class. 
See <?php $filepath = $_GET["filepath"];
echo "<a href='MadGraph5Processes.php?filepath=".$filepath."' target='page'>";?>MadGraph5 Processes</a> 
for a brief overview, and <code>examples/main34.cc</code> for an example 
how to use it. Still at an experimental stage, and only tested for 
Madgraph5_aMC@NLO v2.3.3.</li> 
 
<li>By default the program will now assign the PYTHIA mass for 
massless <i>c</i> and <i>b</i> quarks in Les Houches input, see 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>LesHouches:setQuarkMass</a></code>. 
Usually a logical recoiler is found that can transfer the needed 
four-momentum to the quark. The already existing machinery for 
giving masses to massless leptons has been expanded to use the more 
sophisticated algorithm now in place.</li> 
 
<li>Significant bug found in timelike and spacelike showers, whereby 
the azimuthal anisotropy from gluon polarization in the past has been 
overestimated. This does not affect multijet rates, but can influence 
distributions sensitive to angular correlations, although checks have 
not revealed any appreciable effects. New switches 
<code>TimeShower:phiPolAsymHard</code> and 
<code>SpaceShower:phiPolAsymHard</code> regulate whether the first 
branching after/before the hard process can correlate with the 
hard-process event plane. Thanks to Radek Zlebcik.</li> 
 
<li>Some further new options and minor additions in the machinery for 
hard diffraction. This includes three options for setting the impact 
parameter for the <i>Pomeron p</i> subcollision, a possibility to 
access both it and the impact parameter both for the original 
<i>p p</i>, and options to rescale the Pomeron flux, including one 
that uses the MBR renormalization. Some default values changed, notably 
that now MPI is checked.</li> 
 
<li>A new constructor for the <code>Pythia</code> class takes references 
to a <code>Settings</code> and a <code>ParticleData</code> object as 
inputs. In cases where multiple <code>Pythia</code> copies are created 
this allows the  <code>xmldoc</code> files to be read only once. 
Updated <code>examples/main19.cc</code> illustrates.</li> 
 
<li>New method <code>Particle::daughterListRecursive()</code> 
that uses the <code>daughterList()</code> method to trace 
consecutive generations of decay products.</li> 
 
<li>New <code>Vec4</code> friend method <code>pShift(...)</code> 
to transfer four-momentum between two four-vectors so as to bring them 
to have specified new masses.</li> 
 
<li>New particle method <code>intPol()</code> returns the polarization 
as an integer if the stored double-precision number is very close to 
0, +-1, +-2 or 9, and else -9.</li> 
 
<li>Initialize the random number generator earlier, so a non-default 
seed choice also could benefit early external initialization making 
use of it.</li> 
 
<li>Minor fix in interface to zlib.</li> 
 
<li>Changed default setting in <code>main89mlm.cmnd</code>, to better 
agree with common practice.</li> 
 
<li>Minor improvements and fixes in the weighting facilities for 
initial- and final-state showers.</li> 
 
<li>Minor update in the beam-remnant handling for DIS.</li> 
 
<li>Minor improvements in the handling of resonance mass selection.</li> 
 
<li>The GetDJR function of the <code>JetMatchingMadgraph</code> class 
has been renamed <code>getDJR</code> to adhere to standard naming 
conventions. A pointer in the same class is explicitly nulled.</li> 
 
<li>Bug fix in the selection of masses in resonance decays. 
In rare situations this could give wrong masses for particles. 
Thanks to Are Raklev and Anders Kvellestad.</li> 
 
<li>The <code>StringFlav::combine( int, int bool)</code> method is 
renamed <code>combineId</code> to avoid a potential incorrect 
method overloading. Thanks to James Monk.</li> 
 
<li>Bug fix: copy vertex information when a long-lived particle 
decays to three quarks (typically with baryon number violation), 
whereof two have such a small invariant mass that they collapse to a 
diquark. Thanks to Cristiano Alpigiani.</li> 
 
<li>Bug fix for excited quarks <i>q^*</i> and leptons <i>l^*</i>, 
that if new decay channels were introduced they could incorrectly make 
use of the matrix element expressions for the existing decay modes. 
Thanks to Simone Amoroso.</li> 
 
<li>Bug fix in the kinematics of four or more resonance decay products 
when  kinematics is redone owing to matrix-element corrections. 
Thanks to Simone Amoroso.</li> 
 
<li>Changed off-shell behaviour for squark pair production.</li> 
 
<li>Minor fix for random number start-up in the PowhegBox interface, 
and inserted warning that using LHAPDF5 in both Pythia and PowhegBox 
can be dangerous.</li> 
 
<li>Correct two misspelt endtags for LHEF3 output from Pythia. 
Minor technical changes in the LHEF3 machinery.</li> 
 
<li>Bug fix for information on the pdf value chosen for the hardest MPI, 
which was reported a factor 9/4 too large for an incoming gluon. Does 
not affect the event generation itself.</li> 
 
<li>Correct <code>BeamRemnants:primordialKThard</code> from 2.0 to 1.71 
for ATLAS tune AZ. Thanks to Christian Bauer.</li> 
 
<li>Introduce protection against (close-to-)zero-energy partons in string 
length calculations, and against topologies with extremely small angles 
between two junction legs. Thanks to Jan Fiete Grosse-Oetringhaus.</li> 
 
<li>Check by trial whether a given LHAPDF5 set contains photons or not, 
therey avoiding explicit enumeration of such sets.</li> 
 
<li>When the ARCH environment variable does not have a valid value it 
is set to LINUX in configure. Thanks to Alessandro Degano.</li> 
 
<li>Make the BOOST include directory available when the LHAPDF6 plugin 
is used.</li> 
 
<li>Minor corrections for the LHAPDF6 description in README. 
Thanks to Radek Zlebcik.</li> 
 
<li>Fix harmless name overloading in <code>ProcessContainer.cc</code>.</li> 
 
<li>Updated address for Philip Ilten.</li> 
 
<li>A few bibliography updates.</li> 
 
<li>Year updated to 2016.</li> 
 
</ul> 
</li> 
 
<li>8.212: 23 September 2015 
<ul> 
 
<li>Included new weighting facilities in initial- and final-state showers. 
This allows to consistently enhance rare shower branchings without 
impairing the no-emission probabilities, see 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a>. A new <code>main63.cc</code>, 
together with <code>main63.cmnd</code>, illustrates this feature.</li> 
 
<li>Added new functionality to time- and spacelike showers to streamline 
merging with external shower plugins, see 
<?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>Implement New Showers</a>.</li> 
 
<li>New possibility to access hadronization parameters in each step of 
the hadronization process, and to veto individual hadrons, see 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a>. Thanks to Christian Bierlich. 
</li> 
 
<li>The <?php $filepath = $_GET["filepath"];
echo "<a href='ExternalDecays.php?filepath=".$filepath."' target='page'>";?>EvtGen interface</a> has been 
expanded with the possibility to force a rare decay of some of the 
particle species handled by EvtGen, with an event weight compensating 
for this bias. New status codes 95 or 96 single out particles from a 
forced decay, without or with oscillation. Documentation has been 
expanded.</li> 
 
<li>New option <code>MultipartonInteractions:pTmaxMatch = 3</code> 
introduced to allow better matching for multiparton final states. 
Thanks to Paolo Gunnellini.</li> 
 
<li>Added possibility to write Les Houches Event Files abiding to the 
latest (v. 3) standard,see 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a>. The new 
<code>main64.cc</code> illustrates how to use the LHEF3 writer.</li> 
 
<li>Fixes in the merging machinery, related to the scales associated 
with final-state splittings with an initial-state recoiler. The 
evolution scales and energy sharing values associated with this 
configuration were not calculated correctly, as an incorrect recoiler 
momentum was used, leading to a numerical small deviation from the 
actual result. This issue has now been corrected. Some additional 
matching and merging changes, mostly to make the reclustering in the 
merging numerically more stable, and to make UN2LOPS possible.</li> 
 
<li>Fix bug for the <code>Beams:newLHEFsameInit = on</code> option, 
whereby several LHE files can be read in without new initialization 
for each. When the gzip support was improved in 8.209 the filenames 
were no longer reset properly, and so no files beyond the first one 
were read. Thanks to Roberto Franceschini.</li> 
 
<li>Bug fix in the junction handling for the case of gluinos.</li> 
 
<li>Avoid segmentation fault for some rare junction topologies. 
Thanks to Gabriel Magill.</li> 
 
<li>Minor fix for the path to the PowhegBox plugin, as used in the 
<code>Makefile</code>.</li> 
 
<li>Minor fix for the <code>configure</code> script, for the case 
of building a shared library on Mac OS X. Thanks to Mikhail Kirsanov. 
</li> 
 
<li>Several minor changes in the hard diffractive machinery.</li> 
 
<li>PDG codes 51 - 60 have been reserved for Dark Matter (DM) particles. 
51, 52 and 53 are intended for spin 0, 1/2 and 1 stable DM particles, 
respectively. 54 and 55 are for spin 0 and 1 unstable mediators of 
s-channel DM pair creation or annihilation. The rest can be used freely, 
say in models with several DM particles. Masses should be set 
appropriately. Particles 54 and 55 come equipped with a fictitious mix 
of potential decay modes, to be modified according to the specific model. 
By default all have separate antiparticles, with negative codes. 
These can be removed by setting the antiparticle name to be "void". The 
<code>ParticleData::antiName(...)</code> method has been fixed to accept 
this, just like <code>ParticleData::names(...)</code> already did. </li> 
 
<li>Several processes have been set up with wrong process codes, and are 
now fixed. The corrections are 4003 &rarr; 4203, 4004 &rarr; 4204, 
4005 &rarr; 4205, 961 &rarr; 1061, 3124 &rarr; 3123, 3126 &rarr; 3124, 
3144 &rarr; 3143, and 3146 &rarr; 3144.</li> 
 
<li>For <i>tau</i> decay spin correlations the 
<i>f fbar &rarr; V &rarr; f fbar</i> matrix element, 
<i>V = gamma/Z^0</i> or <i>W^+-</i>, is now only used when the 
incoming fermion pair has no other daughters. Else the simpler 
<i>V &rarr; f fbar</i> matrix element is used.</li> 
 
<li>Minor technical changes, to allow for external <i>pT</i>-unordered 
parton showers. No effect for the internal showers.</li> 
 
<li>The SLHA reading operation has become significantly less verbose: 
where previously identical warning messages could be issued for 
many particle codes, now these codes are collected into a single list 
that is written once.</li> 
 
<li>Bug fix in the input of SLHA decay tables. Previously, the first DECAY 
entry could be ignored under certain conditions, so that the internal SUSY 
machinery was used to calculate decay properties. This was the case, in 
particular, when internal SUSY processes were invoked to generate events. 
</li> 
 
<li>Minor correction in the cross section for <i>f fbar &rarr; G*</i> 
for extra-dimensional processes.</li> 
 
<li><code>main85.cc</code> has been updated so it can handle weighted 
events.</li> 
 
<li><code>main89mlm.cmnd</code> contained misleading instructions which 
are now changed.</li> 
 
<li>Protect against rare but disastrous negative mass-squared by numerical 
roundoff in string fragmentation.</li> 
 
<li>Introduce protection against numerical instability of companion quark 
distribution in <i>x &rarr; 1</i> limit.</li> 
 
<li>Use the <code>PythiaStdlib.h</code> header rather than accessing 
Stdlib directly in <code>SusyLesHouches.h</code>.</li> 
 
<li>Output to <code>cerr</code> replaced by ditto to <code>cout</code>. 
Excepted is the examples main programs, where it serves a pedagogical 
function, while a separate <code>cerr</code> is less relevant for big 
batch runs. (<code>Also FJcore</code> code is excepted).</li> 
 
<li>Give final-state showers used for decays access to the beams; 
not needed for PYTHIA itself but for some plugin showers.</li> 
 
<li>Fix minor typo in <code>README</code>, in instruction for enabling 
64 bit compilation.</li> 
 
<li>Insert two missing <code>std::</code> in 
<code>Pythia8Plugins/GeneratorInputs.h</code>.</li> 
 
<li>Minor update in histogram overflow handling.</li> 
 
<li>Add to the documentation that 
<code>TimeShower:weightGluonToQuark = 1</code> is required for consistent 
aMC@NLO merging.</li> 
 
<li>Clarify documentation on which particles are stable by default.</li> 
 
</ul> 
</li> 
 
<li>8.210: 29 June 2015 
<ul> 
 
<li>Bug fix in CKKW-L merging for LHE files, such that the factorization 
and renormalization scales are set by the <code>SCALUP</code> value if 
the <code>muf2</code> and  <code>mur2</code> LHEF3 attributes have not 
been set, and the user has not set any explicit values. This change 
restores the PYTHIA 8.1 behaviour.</li> 
 
<li>Various technical improvements in the machinery for hard 
diffraction.</li> 
 
<li>Correct quark flavour selection when a string spanned directly 
between two junctions is split up.</li> 
 
<li>Check that SK-I and SK-II colour reconnection machineries only 
are called for event topologies they are set up to handle.</li> 
 
<li>Bug fixes in partial widths of the <i>W'</i> boson. Results are 
correct when the <i>W'</i> is a simply rescaled copy of the <i>W</i>, 
but not for more general couplings. Thanks to Mihail Chizhov. </li> 
 
<li>Minor fix in default location of PDF data files in the constructors. 
No practical consequence since correct non-default values are used.</li> 
 
<li>Tiny fix in the <code>configure</code> script, so that CXX options 
containing an equal sign are parsed correctly.</li> 
 
</ul> 
</li> 
 
<li>8.209: 25 May 2015 
<ul> 
 
<li>The "An Introduction to PYTHIA 8.2", arXiv:1410.3012 [hep-ph], 
has now been published in Comput. Phys.Commun. 191 (2015) 159.</li> 
 
<li>The new QCD-based colour reconnection machinery has been expanded 
with checks whether reconnecting dipoles are causally connected. 
This cannot be defined exactly, and therefore several different 
options are available. Also some other minor changes in this part 
of the code, including updates of parameter default values.</li> 
 
<li>New options for <code>ColourReconnection:flipMode</code> allows 
to sidestep the gluon-move handling but still retain the final 
flip step.</li> 
 
<li>New <code>ColourReconnection:forceResonance</code> switch 
allows an additional colour reconnection step after late resonance 
decays. This is especially relevant for 
<i>H^0 &rarr; W^+ W^- / Z^0 Z^0  &rarr; q_1 qbar_2 q_3 qbar_4</i>, 
since the Higgs is so long-lived that its decay is well separated 
from the rest of the event.</li> 
 
<li>The old SK I and SK II colour reconnection models are now available. 
These are specially aimed at the processes 
<i>e^+ e^- &rarr; W^+ W^- / Z^0 Z^0</i>, but are also relevant for 
<i>H^0 &rarr; W^+ W^- / Z^0 Z^0</i>. </li> 
 
<li>Colour reconnection also made possible, optionally, for the 
<code>forceHadronLevel</code> method.</li> 
 
<li>A new switch <code>POWHEG:QEDveto</code> has been introduced to 
steer the treatment of non-QCD radiation in the POWHEG implementation 
of <code>include/Pythia8Plugins/PowhegHooks.h</code>.</li> 
 
<li>The need to link to the Boost library to read gzipped files 
has been eliminated by including iostream classes wrapping the zlib 
compression library, see new files <code>Streams.h</code> and 
<code>Streams.cc</code>.</li> 
 
<li>Linkage to LHAPDF6 so far has been based on the LHAPDF5 compatibility 
mode, which requires no Boost headers. The new configure option 
<code>--with-lhapdf6-plugin=LHAPDF6.h</code> uses native mode, and 
then requires Boost headers.</li> 
 
<li>The LHAPDF6 PDF values are frozen at the respective <i>x</i> 
and <i>Q^2</i> boundary if the <i>(x, Q^2)</i> pair falls outside 
the fit region.</li> 
 
<li>New functions to extract the fit boundaries, the current 
<i>alpha_s</i> value, and the quark masses from LHADPF6.</li> 
 
<li>New option in the HepMC interface, whereby PYTHIA particles can be 
appended to an existing HepMC event. Thanks to Mikhail Kirsanov.</li> 
 
<li>A new runtime interface to the POWHEGBOX matrix element programs, 
bypassing the need for intermediate LHE files. The new files 
<code>include/Pythia8Plugins/LHAPowheg.h</code> and 
<code>include/Pythia8Plugins/PowhegProcs.h</code> contain the 
LHAup class wrapper used to build the POWHEG plugin libraries and 
the simple class that facilitates loading the POWHEG plugins, 
respectively. The new <code>examples/main33.cc</code> demonstrates 
how to use these plugins, and <code>examples/main33.cmnd</code> contains 
the commands needed for POWHEGBOX to run the example.</li> 
 
<li>When reading in particle data from the SLHA interface, changes done 
by the user takes precedence over the SLHA input ones. To be more 
specific, particle data changes by the <code>Pythia::readString</code> 
and <code>readFile</code> methods are buffered and repeated after 
the SLHA initialization.</li> 
 
<li>New <code>#define PYTHIA_VERSION 8.2xx</code> in <code>Pythia.h</code> 
allows user-code preprecessors to make version-specific choices, and 
allows the <code>Pythia</code> class constructor to check that the 
header-file version number matches those of the source code and the 
XML files. Thanks to Pere Mato.</li> 
 
<li>User-defined semi-internal processes can now be accompanied by 
user-defined phase-space generators, via a second optional argument to 
<code>Pythia::setSigmaPtr(SigmaProcess*, PhaseSpace* = 0)</code>. 
Default is the old behaviour, with PYTHIA selecting the phase-space 
generator itself, based on the process type. Alternatively, the user 
may provide a pointer to an instance of an object inheriting from 
the <code>PhaseSpace</code> class. Sufficiently well-tested and general 
such generators could be communicated to the PYTHIA authors for possible 
inclusion in a future release.</li> 
 
<li>Updated initialization of LHEF files in <code>LesHouches.cc</code> 
to ignore file contents enclosed in comment tags, 
<code>&lt;!-- .. --&gt;</code> and/or the special CDATA statement, 
<code>&lt;![CDATA[ .. ]]&gt;</code>. The latter occurs, e.g., in LHEF 
files produced by CalcHep. (It is generally used to store content that 
contains XML-illegal characters like "&lt;" or "&", such as JavaScript 
source.) Thanks to A. Belyaev and A. Pukhov for help with this update.</li> 
 
<li>New flag <code>LesHouches:matchInOut</code>, by default on, 
to recalculate the energies and longitudinal momenta of the incoming 
particles from the outgoing ones for Les Houches input. Reduces effect 
of numerical inconsistencies in input.</li> 
 
<li>Updated EvtGen interface, to allow a pointer to an FSR engine to be 
passed. Thanks to Torben Ferber.</li> 
 
<li>Spin information for <i>tau</i> leptons now also set up for 
<i>Z'</i> and <i>W'</i> decays.</li> 
 
<li>New method<code>AlphaStrong::setThresholds(...)</code> allows to 
set the charm, bottom and top flavour-threshold masses used for the 
running of <i>alpha_strong</i>.</li> 
 
<li>New option for the <code>Event::list()</code> methods allows to show 
momenta with more decimal digits.</li> 
 
<li>New <code>Particle::isFinalPartonLevel()</code> method to tell whether 
a particle belonged to the final state on the parton level of generation 
or not. New <code>main73.cc</code> example illustrates usefulness.</li> 
 
<li>The <code>[]</code> operator is implemented for the <code>Vec4</code> 
class to return its components by index.</li> 
 
<li>New method <code>string Settings::output(string key)</code> returns 
the value of a variable as a string. In <code>Pythia::readString()</code> 
or <code>readFile()</code> calls this can now be used to print a 
current setting value by the command <code>key = ?</code>.</li> 
 
<li>The role of the <code>pythia.forceTimeShower(...)</code> method is 
better explained in the hadron-level standalone documentation, and 
<code>main21.cc</code> has been extended with an example.</li> 
 
<li>The <code>TimeShower::enhancePTmax()</code> method is made virtual.</li> 
 
<li>The <code>iTopCopyId()</code> and <code>iBotCopyId()</code> methods 
now scan all mothers/daughters, rather than only the first and last, 
when searching for a unique flavour match.</li> 
 
<li>New example <code>main62.cc</code> illustrates how a user hook can 
steer the selection of angles in a resonance decay.</li> 
 
<li>Minor correction in the treatment of the highest multiplicity in 
FxFx jet matching.</li> 
 
<li>Minor correction in the interface to aMC@NLO for dijet and 
photon+jet processes.</li> 
 
<li>Inserted missing endtag that corrupted the <code>Tunes.php</code> 
page. Thanks to Tim Martin.</li> 
 
<li>Minor fix for polarization sign when the <i>tau</i> polarization 
is forced.</li> 
 
<li>The meaning of the <code>HiggsXX:parity</code> options for CP mixing 
has been slightly modified and is better described.</li> 
 
<li>Clarification in the documentation that impact-parameter-enhancement 
factor calculation for two hard processes does not work for the 
<i>x</i>-dependent impact-parameter profile option.</li> 
 
<li>Fixed a factor <i>sqrt(2)</i> error in couplings for 
chargino + squark pair production. Thanks to Michihisa Takeuchi.</li> 
 
<li>Some minor bug fixes in the SUSY code, and some speed optimization 
suggested by Martin White.</li> 
 
<li>Fixed incorrect indexing in the SLHA interface.</li> 
 
<li>Removed misleading flavour setup for unused Pomeron.</li> 
 
<li>Fixed possible unwarranted destruction of pointer to external timelike 
shower by the <code>Pythia</code> destructor.</li> 
 
<li>Initialized pointers to NULL in the <code>Info</code> class, 
to avoid some problems. Thanks to Keno Fischer.</li> 
 
<li>Minor <code>Makefile</code> improvements for LHAPDF linking.</li> 
 
<li>Minor improvement to the <code>main89.cc</code> example program.</li> 
 
<li>Mildly modified warning/error messages when junction splitting 
fails.</li> 
 
<li>Minor fixes in <code>LHAFortran.h</code>, which is also moved to 
<code>include/Pythia8Plugins</code> to better reflect its peripheral 
role.</li> 
 
<li>New machinery for hard diffraction now in place. Still being 
debugged and tested, so not yet ready for public usage.</li> 
 
</ul> 
</li> 
 
<li>8.205: 23 January 2015 
<ul> 
 
<li>Unfortunate tiny typo made the new A14 tunes inaccessible. 
Now fixed.</li> 
 
<li>Resolved an inconsistency between MLM and FxFx merging that was 
introduced in 8.201 upon answering a user request. The effect should 
be minimal. Thanks to Josh Bendavid for bringing this to our attention. 
</li> 
 
<li>Amended the automated reading of the (optional) 
"TimeShower:nPartonsInBorn" setting from Les Houches Event files 
produced with aMC@NLO. Heavy coloured objects are now handled correctly 
when the information in the Les Houches event does not include heavy 
partons of the lowest-multiplicity state into its counting (as is 
conventional). Thanks to Josh Bendavid for pointing this out.</li> 
 
</ul> 
</li> 
 
<li>8.204: 22 January 2015 
<ul> 
 
<li>The <code>examples</code> directory has been moved back from 
the <code>share/Pythia8/</code> subdirectory to the main directory, 
as was the case in PYTHIA 8.1, to make it more visible to newcomers. 
The optional <code>make install</code> step will create a copy of 
<code>examples</code> in <code>share/Pythia8/</code>. The rarely 
used <code>examples/outref</code> subdirectory is moved to 
<code>share/Pythia8/outref</code>.</li> 
 
<li>Fifteen new <?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>tunes</a> have been added, 
the MonashStar tune from CMS and fourteen A14 tunes from ATLAS 
[<a href="Bibliography.php" target="page">ATL14a</a>]. The latter correspond to central tunes for 
four different PDF sets and ten variations in five (approximate) 
eigenvector directions. Furthermore, now the chosen 
<code>Tune:pp</code> implies the <code>Tune:ee</code> value 
to which it is related, and thus the latter need not be set 
separately.</li> 
 
<li>Default settings values have been updated to agree with the 
Monash 2013 tune. Thus typically the list of changed settings is 
significantly reduced. Thanks to Mikhail Kirsanov for suggestion.</li> 
 
<li>The compositeness section has been expanded with six further 
processes, describing the pair production of excited leptons or 
neutrinos. Three-body contact-interaction decay modes of these 
excited states have been introduced. A bug has been fixed that 
gave the wrong helicity in decays for excited quarks, leptons and 
neutrinos. Further, the <i>gamma^*/Z^0/Z'^0</i> can decay to a pair 
of excited fermions provided that the channels are added to the 
list of allowed ones. Based on code provided by Olga Igonkina.</li> 
 
<li>A new interface to the EvtGen decay package, primarily intended 
for bottom and charm decays, has been implemented. It is available 
in <code>include/Pythia8Plugins/EvtGen.h</code> and an example how 
to set it up is found in <code>main48</code>.</li> 
 
<li>Several minor changes in the <code>examples</code> subdirectory. 
This includes the <code>README</code>, <code>Makefile</code> an 
<code>runmains</code> files. The <code>main61</code> example program 
has been removed, since now LHAPDF can be loaded dynamically for 
<code>main42</code>, so that the two become equivalent. Further 
<code>main62</code> has been renamed <code>main43</code>, to gather 
HepMC-related examples, and <code>testlhef3.lhe</code> has been 
renamed <code>wbj_lhef3.lhe</code>. 
</li> 
 
<li>A new example <code>main30</code> how to create a tailormade 
copy of the ordinary event record, here with a history tracing 
of the hard process closer to the PYTHIA 6 conventions.</li> 
 
<li>The <?php $filepath = $_GET["filepath"];
echo "<a href='ProMCFiles.php?filepath=".$filepath."' target='page'>";?>ProMC</a> input-output file format 
is now implemented among the libraries that can be 
<code>configure</code>d to run with PYTHIA. An examples is provided 
in <code>main46.cc</code>. Thanks to Sergei Chekanov.</li> 
 
<li>Change in the setup of final-state-shower colour dipoles for the 
non-default case of no interleaving, whereby it becomes less likely 
to pick a colourless final-state particle as recoiler. New option 
<code>TimeShower:allowMPIdipole</code> gives more flexibility. 
Thanks to Mihoko Nojiri and Bryan Webber.</li> 
 
<li>New options 3 and 4 for <code>TimeShower:pTdampMatch</code> 
and <code>SpaceShower:pTdampMatch</code>, with new default 3 for the 
latter. The main effect is that, by default, <i>t tbar</i> 
production (as a <i>2 &rarr; 2</i> process) obtains damped 
radiation above the process scale. Thanks to Andy Buckley 
for suggestion.</li> 
 
<li>Initialization will now abort if a mode has been chosen with a 
non-allowed value. This applies to those modes that have been 
defined with the <code>modepick</code> or <code>modefix</code> 
labels in the <code>xmldoc/*.xml</code> files and, for the former, 
where maximal and minimal values have been specified. The former 
label is used to represent a discrete set of options, and so any 
value outside the allowed range is just plain wrong. Thanks to James 
Monk for suggestion.</li> 
 
<li>Added <code>getChannels()</code> in <code>SusyResonanceWidths</code> 
to dynamically create the decay table for SUSY particles and thereby 
to remove duplication in the XML file.</li> 
 
<li>New file <code>SusyWidthfunctions</code>. Added new class 
<code>WidthFunction</code> to handle calculation of three- and four-body 
decay widths.</li> 
 
<li>In <code>main24.cc</code> the example spectrum file format has been 
updated from SLHA1 to SLHA2, obtained from SoftSusy 3.5.1. The two new 
<code>slha1-example.spc</code> and <code>slha2-example.spc</code> files 
replace the older <code>cmssm.spc</code>, <code>snowmass2.spc</code> 
and <code>softsusy.spc</code> ones. </li> 
 
<li>In <code>LesHouches.cc</code> the read-in of LHEF headers 
containing tags that open and close on the same line (e.g., a single 
line containing <code>&lt;tag&gt;blabla&lt;/tag&gt;</code>) has been 
enabled. This could previously lead to improper initialization and 
crashes. Also implemented a check for forgotten close-tag statements, 
with a warning issued to the user. Thanks to Alexander Belyaev and 
Alexander Pukhov.</li> 
 
<li>Bug fix in the conversion from the <code>xml</code> settings files 
to their <code>php</code> radio-button equivalents, whereby the text 
describing some options was not set properly. (Typically the text of 
the previous option was repeated.) Thanks to Radek Zlebcik.</li> 
 
<li>Minor fixes in the LHEF version 3 reader. Introduce a new 
matching writer of LHEF version 1 or 3 files.</li> 
 
<li>Introduction of a new mode <code>LesHouches:setLeptonMass</code>, 
such that by default final-state charged leptons acquire sensible 
masses, even when the matrix-element calculations have been 
performed with massless leptons. This and other energy-momentum 
adjustments, e.g. for limited-precision storage, are now located in 
<code>ProcessContainer::constructProcess</code>.</li> 
 
<li><code>http://home.thep.lu.se/Pythia</code> has been introduced 
as a simpler but (hopefully) equivalent address to 
<code>http://home.thep.lu.se/~torbjorn/Pythia.html</code>, and 
various documentation has been updated accordingly. Thanks to 
Leif L&ouml;nnblad.</li> 
 
<li>New file <code>include/Pythia8Plugins/execinfo.h</code> contains 
trivial copies of three backtrace methods needed to be able to compile 
PYTHIA under Cygwin. The <code>README</code> file and the worksheet 
updated with brief information on three (non-supported) ways of working 
with PYTHIA under Windows. Thanks to and Theo Hughes and Gordon Watts. 
</li> 
 
<li>Bug fix in check for colour sextets and transfer of such colour 
information. Thanks to Alexander Belyaev and Alexander Pukhov.</li> 
 
<li>Improved handling of stray characters in the SUSY Les Houches 
code. The check on the consistency of decay tables has been removed. 
Improved warning/error printing in the SLHA interface.</li> 
 
<li>Bug fix in new beam remnant model, so that it basically 
operates like the old one for <i>e^+e^-</i> annihilation.</li> 
 
<li>Two bug fixes in the new colour reconnection model, one for 
diquarks at the ends of junction strings, and another to check that 
coloured resonances are processes with early resonance decays option.</li> 
 
<li>Bug fix for multiple <code>Pythia::init()</code> calls, where 
beam contents were not properly reset. Thanks to Josh Bendavid.</li> 
 
<li>Bug fix such that the valence content of a <i>pi^0</i>, 
<i>K^0_S</i>, <i>K^0_L</i> and Pomeron is reselected for each 
new event. Thanks to Radek Zlebcik.</li> 
 
<li>Fix typo in constants of the <i>tau &rarr; 3 pi</i> current 
for the amplitudes of the <i>rho</i>, <i>rho(1450)</i>, and 
<i>f2</i>. Thanks to Ian Nugent.</li> 
 
<li>Small bug fixes for string and ministring fragmentation, for the 
case when a low-mass (order 2 GeV) system contains at least three 
partons, which fail to define a unique direction for the final 
string region.</li> 
 
<li>New parameter <code>BeamRemnants:reducedKTatHighY</code> introduced 
to reduce technical problems with low-mass MPIs produced at high 
rapidities when primordial <i>kT</i> is introduced.</li> 
 
<li>Small bug fix in the global-recoil option for timelike showers.</li> 
 
<li>Update year to 2015, remove tabs and superfluous blanks, break long 
lines where meaningful, and some further minor changes.</li> 
 
</ul> 
</li> 
 
<li>8.201: 14 October 2014 
<ul> 
 
<li>The <i>Introduction to PYTHIA 8.2</i> has now been assigned 
the arXiv:1410.3012 [hep-ph] identifier, which has been introduced 
in code and text.</li> 
 
<li>The <code>enable-shared</code> by mistake was not listed 
among allowed configure options.</li> 
 
<li>Corrected a few tiny documentation typos.</li> 
 
</ul> 
</li> 
 
<li>8.200: 11 October 2014 
<ul> 
 
<li>A new <code>share/Pythia8</code> directory collects all 
documentation and example code. The <code>examples</code>, 
<code>htmldoc</code>, <code>phpdoc</code> and <code>xmldoc</code> 
directories have been moved here. The main-directory files 
<code>AUTHORS</code>, <code>COPYING</code>, <code>GUIDELINE</code> 
and <code>README</code> are also copied here during installation. 
</li> 
 
<li>A new <code>share/Pythia8/pdfdoc</code> directory collects pdf 
documents that are linked from the <code>htmldoc</code> and 
<code>phpdoc</code> directories. Over time it will  provide more 
in-depth descriptions of various physics aspects than offered in 
the html/php-formatted documentation. In addition to the official 
main publication and the worksheet, currently notes on LO vs. NLO 
PDFs and on the <i>g &rarr; q qbar</i> branching kernel are 
included.</li> 
 
<li>A new <code>include/Pythia8Plugins</code> directory collects 
code that does not form part of the core PYTHIA functionality but 
still has a general usefulness. Code in this directory will not be 
compiled as part of the Pythia library, but can be linked where needed. 
This new directory contains 
<ul> 
<li>the jet matching classes in <code>CombineMatchingInput.h</code>, 
<code>GeneratorInput.h</code> and <code>JetMatching.h</code>, moved 
from the <code>examples</code> directory;</li> 
<li>the <code>PowhegHooks</code> user hook, to veto shower emissions 
above the POWHEG scale, formerly found in <code>examples/main31.cc</code>; 
</li> 
<li>the <code>Pythia8ToHepMC</code> interface for output of PYTHIA events 
into the HepMC format, combining the code previously in 
<code>include/Pythia8ToHepMC.h</code> and 
<code>pythia8tohepmc/Pythia8ToHepMC.cc</code> into a new 
<code>HepMC2.h</code> file;</li> 
<li>the <code>FastJet3.h</code> interface of PYTHIA particles to the 
FastJet 3 library of jet finders, formerly found in 
<code>include/FastJet3.h</code>; and</li> 
<li>the <code>LHAPDF5.h</code> and <code>LHAPDF6.h</code> files for 
interfaces to the LHAPDF library (see further below).</li> 
</ul></li> 
 
<li>The configure/make structure has been considerably rewritten. 
Now all external libraries to be linked are specified in the 
main-directory <code>configure</code> step, along with other options, 
so there is no longer an <code>examples/configure</code>. The 
<code>make</code> step will, as before, compile and install libraries 
inside the current directory, such that the main programs in the 
<code>examples</code> directory can be run. One small difference is that 
also the archive libraries are installed in <code>lib</code> and not in 
<code>lib/archive</code>. 
<br/>A new optional <code>make install</code> step allows you to copy 
files to more convenient locations. The default option, with no directories 
specified in the <code>configure</code> step, requires you to have 
superuser privileges. Then files will be copied to standard locations 
as follows: 
<table border="0"> 
  <tr> <td>lib/</td> <td>&rarr;&nbsp;</td> <td>/usr/lib/</td> </tr> 
  <tr> <td>include/</td> <td>&rarr;&nbsp;</td> <td>/usr/include/</td> </tr> 
  <tr> <td>share/</td> <td>&rarr;&nbsp;</td> <td>/usr/share/</td> </tr> 
  <tr> <td>pythia-config</td> <td>&rarr;&nbsp;</td> <td>/usr/bin/</td> </tr> 
</table> 
</li> 
 
<li> 
The <code>pythia8-config.in</code> script has been replaced by a new 
<code>bin/pythia8-config</code> script. See the README file for details. 
The <code>make install</code> step by default will put a copy of it in 
<code>/usr/bin</code>. 
</li> 
 
<li>The interface to LHAPDF is now dynamically loaded when requested, 
and can be either to version 5 or 6 of the library. The dummy code 
previously in <code>lhapdfdummy/LHAPDFDummy.cc</code>, to be linked 
when LHAPDF is not, is no longer required. The two new files 
<code>LHAPDF5.h</code> and <code>LHAPDF6.h</code> in the 
<code>include/Pythia8Plugins</code> directory contain the necessary 
interface code. The selection of PDF sets, notably for the proton, 
has been extended to simplify mixing of internal and external PDF sets, 
and it is now possible to specify different PDFs for the two incoming 
protons at the LHC, see the <?php $filepath = $_GET["filepath"];
echo "<a href='PDFSelection.php?filepath=".$filepath."' target='page'>";?>PDF Selection</a> 
description.</li> 
 
<li>The new <code>LHEF3.h</code> file contains a generic interface for 
reading Les Houches Event Files of versions 1.0, 2.0 and 3.0. This 
allows more information to be read and studied by the author. Currently 
PYTHIA itself makes little use of the information beyond the one in 1.0, 
but it is available among the 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Event Information</a>. 
Examples are found in <code>main37.cc</code> and <code>main38.cc</code>. 
</li> 
 
<li>The new <code>Beams:strictLHEFscale</code> switch can be used to 
restrict parton showers in resonance decays to be below the input 
Les Houches scale, not only the hard process itself. 
The new <code>Beams:setProductionScalesFromLHEF</code> switch can be used 
to restrict the emission off each separate parton to be below its specific 
scale.</li> 
 
<li>The <code>rootexamples</code> directory has been removed, and the 
two programs <code>examples/main91</code> and <code>examples/main92</code> 
now illustrate how ROOT can be used in conjunction with PYTHIA.</li> 
 
<li>The executable built from <code>examples/mainxx.cc</code> is now 
named <code>examples/mainxx</code>, while previously it was named 
<code>examples/mainxx.exe</code>.</li> 
 
<li>The rudimentary support for compilation on Windows platforms, 
present in PYTHIA 8.1, has not yet been updated for 8.2 and so is omitted. 
Also the README.HepMC file is omitted for now.</li> 
 
<li>The ProMC interface is broken, and has been removed for now.</li> 
 
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
original input string.</li> 
 
<li>The default tune has been changed from 4C to Monash 2013, meaning 
<code>Tune:ee = 7</code> and <code>Tune:pp = 14</code>. The old 4C 
tune that was default in 8.1 can be recovered with 
<code>Tune:ee = 3</code> and <code>Tune:pp = 5</code>. 
Also most other older tunes are based on <code>Tune:ee = 3</code>. 
</li> 
 
<li>Two new CMS underlying-event tunes [<a href="Bibliography.php" target="page">CMS14</a>] and the ATLAS 
AZ tune [<a href="Bibliography.php" target="page">ATL14</a>] have been added as options.</li> 
 
<li>The default handling of the <i>g &rarr; q qbar</i> splitting kernel 
has been changed, affecting in particular heavy-flavour production. 
<code>TimeShower:weightGluonToQuark</code> has been changed from 1 to 4 
to do this. All old tunes are with the 1 value but, since the tunes are 
not probing the detailed <i>g &rarr; q qbar</i> behaviour, this is 
not set as part of the tune options.</li> 
 
<li>Christine O. Rasmussen joins as new PYTHIA collaboration member.</li> 
 
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
new settings have been introduced, notably 
<code>ColourReconnection:mode</code> to switch among the three models. 
</li> 
 
<li>A new <code>include/Pythia8Plugins/ColourReconnectionHooks.h</code> 
makes available an even larger selection of toy colour reconnection 
models, via user hooks. Some of them are only intended for top decays, 
for top mass uncertainty studies, whereas others can be used more 
generally. The <code>examples/main29.cc</code> program illustrates how 
the different options should be set up.</li> 
 
<li>Several new features and improvements in the matching/merging 
machinery. Notably the aMC@NLO matching scheme has been implemented, 
see the <?php $filepath = $_GET["filepath"];
echo "<a href='aMCatNLOMatching.php?filepath=".$filepath."' target='page'>";?>aMC@NLO Matching</a> 
description. To this end the global-recoil option of timelike showers 
has been improved, and security checks have been introduced for 
inaccurate LHEF input. A new <code>main89.cc</code> example has been 
introduced, where different <code>.cmnd</code> files show how to set 
up either CKKW-L, FxFx, MLM, UMEPS or UNLOPS merging.</li> 
 
<li>Improved capability for the <code>LHAup</code> Les Houches interface 
to read SLHA information embedded in the input file or stream.</li> 
 
<li>The <code>Makefile</code>s have been updated to take into account 
the changed structure of the HepMC interface.</li> 
 
<li>The <i>Z'</i> production process has been updated to optionally 
allow decay to a fourth generation of fermions, with universal or 
non-universal couplings.</li> 
 
<li>Introduction of a new Higgs CP-mixing parametrization via a mixing 
angle <i>phi</i> as described in <?php $filepath = $_GET["filepath"];
echo "<a href='HiggsProcesses.php?filepath=".$filepath."' target='page'>";?>Higgs 
Processes</a>. The choice of the Higgs CP-mixing parametrization 
now also affects the distributions of the <i>tau</i> decay products 
from the processes <i>H^0 &rarr; tau^+ tau^-</i>. 
 
<li>Bug fix in <i>H^0 &rarr; W^+ W^- &rarr; 4 f</i> matrix element 
for mixed CP-state case.</li> 
 
<li>Various improvements and finer grain control for the determination 
of <i>tau</i> decay correlations and <i>tau</i> polarizations. By 
default the decays of <i>tau</i> pairs from known resonance decays 
in Les Houches input are now correlated. 
The <code>ParticleDecays:sophisticatedTau</code> mode 
in <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDecays.php?filepath=".$filepath."' target='page'>";?>Particle Decays</a> has been renamed 
<code>TauDecays:mode</code>, as well as all <i>tau</i>-related 
<code>ParticleDecay</code> options, with two new options of 
using only the internal machinery to determine correlations and 
polarizations, and only using the provided SPINUP digit from Les 
Houches input. The option <code>TauDecays:externalMode</code> has been 
introduced to control the interpretation of the SPINUP digit. 
</li> 
 
<li>For Les Houches Event input the energy of a particle is recalculated 
from its three-momentum and mass, in order to limit mismatches from 
limited numerical precision in the input values.</li> 
 
<li>Bug fix in the two-loop running <i>alpha_s</i>, for the matching 
to six flavours at the top mass.</li> 
 
<li>Eliminate harmless compiler warnings for <code>FJcore</code>.</li> 
 
<li>Updated Introduction (= the official 8.2 article) and Worksheet.</li> 
 
</ul> 
</li> 
 
</ul> 
 
</body>
</html>
 
<!-- Copyright (C) 2016 Torbjorn Sjostrand --> 
