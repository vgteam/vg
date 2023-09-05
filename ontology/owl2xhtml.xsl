<!--
(C) 2016 SIB Swiss Institute of Bioinformatics, http://www.isb-sib.ch
(C) 2014 UniProt consortium, http://www.uniprot.org
(C) 2008 by Andreas Radinger and Martin Hepp, Bundeswehr University Munich, http://www.unibw.de/ebusiness/

This file is part of owl2xhtml, a stylesheet for client-side rendering OWL ontologies in the form of XHTML documentation. For more information, see http://www.ebusiness-unibw.org/projects/owl2xhtml/.
It has been adapted for use on uniprot.org.    

Acknowledgements: The stylesheet re-uses, with kind permission, code-snippets from the RDFS/OWL presentation stylesheet (1) 

by Masahide Kanzaki, and from the OWL2HTML stylesheet (2), by Li Ding. We are very grateful for this kind support.

(1) http://www.kanzaki.com/ns/ns-schema.xsl
(2) available at http://daml.umbc.edu/ontologies/webofbelief/xslt/owl2html.xsl, 


     owl2xhtml is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License (LPGL)
     as published by the Free Software Foundation, either version 3 of 
     the License, or (at your option) any later version.

     owl2xhtml is distributed in the hope that it will be useful, 
     but WITHOUT ANY WARRANTY; without even the implied warranty of 
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
     GNU General Public License for more details. 

     You should have received a copy of the GNU General Public License 
     along with owl2xhtml.  If not, see <http://www.gnu.org/licenses/>. 

-->



<!DOCTYPE xsl:stylesheet [
  <!ENTITY rdf		"http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <!ENTITY foaf "http://xmlns.com/foaf/0.1/" >
  <!ENTITY rdfs		"http://www.w3.org/2000/01/rdf-schema#">
  <!ENTITY dc		"http://purl.org/dc/elements/1.1/" >
  <!ENTITY skos		"http://www.w3.org/2004/02/skos/core#" >  
  <!ENTITY xsd		"http://www.w3.org/2001/XMLSchema#">
  <!ENTITY owl		"http://www.w3.org/2002/07/owl#">
  <!ENTITY core		"http://purl.uniprot.org/core/">
  <!ENTITY vg		"http://biohackathon.org/resource/vg#">
  <!ENTITY spin         "http://spinrdf.org/spin#">
]>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:taxo="http://purl.org/rss/1.0/modules/taxonomy/"
	xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:owl="http://www.w3.org/2002/07/owl#" xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
	xmlns:dc="http://purl.org/dc/elements/1.1/"
>
	<xsl:output method="html" version="1.0" encoding="UTF-8" indent="yes" />

  <!-- global variables and functions -->
	<xsl:variable name="class-name" select="//owl:Ontology/rdfs:label" />
	<xsl:variable name="nodeset-ontology"
		select=".//*[
		  rdf:type/@rdf:resource='&owl;Ontology'  or (local-name()='Ontology' and namespace-uri()='&owl;')
  ]" />
	<xsl:variable name="nodeset-class"
		select=".//*[     
		   (rdf:type/@rdf:resource='&owl;Class'  or (local-name()='Class' and namespace-uri()='&owl;' and @rdf:about) or (owl:Class and @rdf:about))
  ]" />
	<xsl:variable name="nodeset-property"
		select=".//*[  
		(   (local-name()='Property' and namespace-uri()='&rdf;')
		or (local-name()='ConstraintProperty' and namespace-uri()='&rdfs;')
		or (local-name()='DatatypeProperty' and namespace-uri()='&owl;')
		or (local-name()='ObjectProperty' and namespace-uri()='&owl;') 
		or (rdf:type/@rdf:resource='&rdfs;ConstraintProperty')
		or (rdf:type/@rdf:resource='&owl;DatatypeProperty') 
		or (rdf:type/@rdf:resource='&owl;ObjectProperty')
		or (rdf:type/@rdf:resource='&rdf;Property') 
		)]" />
	<xsl:variable name="nodeset-property-object" select=".//*[  
		((local-name()='ObjectProperty' and namespace-uri()='&owl;') or (rdf:type/@rdf:resource='&owl;ObjectProperty')) 
  ]" />
	<xsl:variable name="nodeset-property-datatype" select=".//*[  
		((local-name()='DatatypeProperty' and namespace-uri()='&owl;') or (rdf:type/@rdf:resource='&owl;DatatypeProperty'))
  ]" />
	<xsl:variable name="nodeset-individual"
		select="(.//*[
			 (//rdf:type[@rdf:resource='&owl;NamedIndividual'] or (//rdf:type[@rdf:resource='&owl;Thing'] and (namespace-uri()='&core;')))
		  and not (namespace-uri()='&rdfs;')
		  and not (local-name()='RDF')
		  and @rdf:about
  ])" />
	<xsl:template match="/">
	<!-- Put IE into quirks mode -->
		<html>
			<head>
				<title>
					<xsl:value-of select="//owl:Ontology/rdfs:label" />
				</title>
				<style type="text/css">
                                  .Class {padding:0.4em; width:99%; border-width: 1px; border-style: dashed; position:relative;
					margin-top: 1.5em; border-style:solid; background:#e0e0e0}
                                  .cp-type {font-weight: normal; font-size:90%; color:
					maroon;} .nice-content td { width:50% }
				  #header {
                                        padding-top: 2px;
                                        padding-left: 1px;
                                        padding-right: 1px;
					background-color: #333333;
					
					color: #ffffff;
				  }
				  #header a {
                                        color: #ffffff;
                                  }
				  #header #logo {
                                     padding-left: 20px;
				  }
				  #header #trail {
                                     text-align: left;
				  }
				  #header #menu {
                                     text-align: right;
                                     padding-right: 20px;
				  }
				  body {
                                    margin:0;
                                    padding:0;
                                    font-family: trebuchet ms,sans-serif!important;
				  }
				  table {
                                    width:100%;
                                    border-collapse: collapse;
				  }
				  #search {
                                      margin-top:1em;
                                      margin-left: 0.5em;
                                      padding-left: 1px;
                                      padding-right: 1px;
                                      background-color: #F5F5F5;
				  }
				  h2, h3 {
                                    border-bottom: 2px solid #ACBCD0;
                                    color: #426490;
                                    font-size: medium;
                                    font-weight: bold;
                                    margin-bottom: 0.5em;
                                    margin-top: 1em;
                                    padding: 0.25em 0;
                                  }
                                  .subsection {
                                    margin-top: 0.5em;
                                  }
                                  th {
                                     text-align: left;
                                  }
                                  .nice {
                                    color: #000000;
                                  }
                                  .nice th, .nice td {
                                      padding: 0.2em 0.3em;
                                      vertical-align: top;
                                  }
                                  .nice tr td {
                                    border-left: 1px solid #DEDEDE;
                                  }
                                  .nice tr td:first-child  {
                                    border-left: 0px solid #DEDEDE;
                                  }
                                  .nice tr {
                                    border-top: 1px solid #DEDEDE;
                                  }
                                  .nice tr:first-child {
                                    border-top: 0px solid #DEDEDE;
                                  }
                                  #content {
                                    margin-left:0.5em;
                                    margin-right:0.5em;
                                  }
                                  .indent {
                                    margin-left:2em;
                                  }
                                  #footer {
                                    margin-top:1em;
                                    padding-left:1px;
                                    padding-right:1px;
                                    background-color: #F5F5F5;
                                    margin-bottom:0em;
                                  }
                                  #footer ul {
                                     list-style: none outside none;
                                  }
				</style>
			</head>
			<body>
				<xsl:apply-templates />
				<div id="footer">
                                    <p>We would like to especially thank the organisers and funders of the <a href="http://www.biohackaton.org">BioHackathon</a> series of meetings for hosting the original discussions leading to VG RDF</p> 
                                    <p>VG RDF is addittionaly supported by:</p>
                                    <ul>
                                      <li><a href="http://biosciencedbc.jp/en/">National Bioscience Database Center (NBDC)</a> of <a href="http://www.jst.go.jp/">Japan Science and Technology Agency (JST)</a></li>
                                      <li><a href="http://dbcls.rois.ac.jp/en/">Database Center for Life Science (DBCLS)</a></li>
                                      <li><a href="http://www.rois.ac.jp/">Research Organization of Information and Systems (ROIS)</a></li>
                                      <li><a href="http://www.isb-sib.ch/groups/geneva/sp-xenarios.html">Swiss-Prot part of the SIB Swiss Institute of Bioinformatics</a> and supported by the <a href="http://www.sbfi.admin.ch/?lang=en">Swiss Federal Government through the The State Secretariat for Education, Research and Innovation SERI</a>.</li>
                                      <li><a href="http://www.nih.gov/">NIH under R24OD011883</a> </li>
                                      <li><a href="http://www.ddbj.nig.ac.jp">DNA Databank of Japan (DDBJ)</a></li>
                                      <li><a href="http://www.sanger.ac.jp">Sanger centre</a></li>
                                      <li><a href="https://www.u-tokyo.ac.jp/en/index.html">The University of Tokyo</a></li>
                                      <li><a href="https://genomics.soe.ucsc.edu/haussler">University of California, Santa Cruz, David Haussler lab</a></li>
                                    </ul>
				</div>
			</body>
		</html>
	</xsl:template>

	<!-- Match RDF data -->
	<xsl:template match="rdf:RDF">
		<table id="header">
			<tr>
				<td id="logo">
					<a href="/" accesskey="1">
						<img src="http://www.biohackathon.org/bh-favicon.png" alt="" width="37px" />
					</a>
				</td>
				<td id="trail">
					VG RDF: RDF properies and classes used to describe genomic variant graphs
				</td>
				<td id="menu">
					<a href="http://biohackathon.org/resource/vg.ttl">Turtle</a>
					&#183; <!-- middot -->
					<a href="http://biohackathon.org/resource/vg.rdf">RDF/XML</a>
                                        &#183; <!-- middot -->
					<a href="https://gitter.im/ekg/vg" accesskey="9">Contact/Help</a>
					&#183; <!-- middot -->
					<a href="https://github.com/ekg/vg" rel="help">Documentation</a>
				</td>
			</tr>
		</table> <!-- end of header -->
                 <table id="search">
                  <tr>
                    <td style="width:20em;">
                            <a href="#OntologyDescription">Ontology description</a>
                    </td>
                </tr>
                <tr>
                    <td>
                            <a href="#Classes">
                                    Classes (
                                    <xsl:value-of select="count($nodeset-class)" />
                                    )
                            </a>
                    </td>
                    <td>
                            <select onChange="window.location.hash = document.getElementById('naviClass').value"
                                            id="naviClass"
                                    >
                                    <xsl:apply-templates select="$nodeset-class" mode="menu">
                                            <xsl:sort select="@rdf:ID" />
                                    </xsl:apply-templates>
                            </select>
                    </td>
            </tr>
            <tr>
                    <td>
                            <a href="#Properties">
                                    Properties (
                                    <xsl:value-of select="count($nodeset-property)" />
                                    )
                            </a>
                    </td>
            </tr>
            <tr>
                    <td>
                            <a style="objectPropertySearch" href="#ObjectProperties">
                                    Object properties (
                                    <xsl:value-of select="count($nodeset-property-object)" />
                                    )
                            </a>
                    </td>
                    <td>
                            <select onChange="window.location.hash = document.getElementById('navi3').value"
                                    id="navi3"
                            >
                                    <xsl:apply-templates select="$nodeset-property-object" mode="menu" />
                            </select>
                    </td>
            </tr>
            <tr>
                    <td>
                            <a style="datatypePropertySearch" href="#DatatypeProperties">
                                    Datatype properties (
                                    <xsl:value-of select="count($nodeset-property-datatype)" />
                                    )
                            </a>
                    </td>
                    <td>
                            <select onChange="window.location.hash = window.location.hash = document.getElementById('navi4').value"
                                            id="navi4"
                                    >
                                            <xsl:apply-templates select="$nodeset-property-datatype" mode="menu">
                                                    <xsl:sort select="@rdf:ID" />
                                            </xsl:apply-templates>
                            </select>
                    </td>
            </tr>
            <tr>
                    <td>
                            <a href="#InstanceData">
                                    Instance data (
                                    <xsl:value-of select="count($nodeset-individual)" />
                                    )
                            </a>
                    </td>
                    <td>
                            <select onChange="window.location.hash = window.location.hash = document.getElementById('navi5').value"
                                            id="navi5"
                                    >
                                            <xsl:apply-templates select="$nodeset-individual" mode="menu">
                                                    <xsl:sort select="@rdf:ID" />
                                            </xsl:apply-templates>
                                    </select>
                    </td>
            </tr>
    </table>
		<div id="content">
			<div id="sections">
				<div class="nice" id="section_OntologyDescription" style="position: relative;">
					<h2>
						Ontology description
					</h2>
					<div class="nice-content" id="OntologyDescription">
						<xsl:if test="count($nodeset-ontology)>0">
							<xsl:apply-templates select="$nodeset-ontology" mode="details_ontology">
							</xsl:apply-templates>
						</xsl:if>
					</div>
				</div>
				<div class="nice" id="section_Classes" style="position: relative;">
					<h2>
						Classes
					</h2>
					<div class="nice-content" id="Classes">
						<xsl:if test="count($nodeset-class)>0">
							<xsl:apply-templates select="$nodeset-class" mode="details">
								<xsl:sort select="@rdf:ID" />
							</xsl:apply-templates>
						</xsl:if>
					</div>
				</div>
				<div class="nice" id="section_Properties" style="position: relative;">
					<h2>
						Properties
					</h2>
					<div class="nice-content" id="Properties">
						<h3 id="ObjectProperties">Object properties</h3>
						<xsl:if test="count($nodeset-property-object)>0">
							<xsl:apply-templates select="$nodeset-property-object" mode="details">
								<xsl:sort select="@rdf:ID" />
							</xsl:apply-templates>
						</xsl:if>
						<h3 id="DatatypeProperties">Datatype properties</h3>
						<xsl:if test="count($nodeset-property-datatype)>0">
							<xsl:apply-templates select="$nodeset-property-datatype" mode="details">
								<xsl:sort select="@rdf:ID" />
							</xsl:apply-templates>
						</xsl:if>
					</div>
				</div>
				<div class="nice" id="section_InstanceData" style="position: relative;">			
					<div class="nice-content" id="InstanceData">
						<xsl:if test="count($nodeset-individual)>0">
							<xsl:apply-templates select="$nodeset-individual" mode="details">
								<xsl:sort select="@rdf:ID" />
							</xsl:apply-templates>
						</xsl:if>
					</div>
				</div>
			</div> <!-- end of sections -->
			<script type="text/javascript"> &lt;!-- UniProt.visibility.initSections('owl'); 
			</script>
		</div> <!-- end of container -->
	</xsl:template>
	<xsl:template match="*" mode="menu">
		
			<xsl:choose>
				<xsl:when test="@rdf:about">
					<option value="{substring-after(@rdf:about, '#')}">
						<xsl:call-template name="prettyUrl">
							<xsl:with-param name="name" select="@rdf:about"/>
						</xsl:call-template>
					</option>
				</xsl:when>
				<xsl:otherwise>
					<option>
						<xsl:value-of select="local-name()" />
					</option>
				</xsl:otherwise>
			</xsl:choose>
		
	</xsl:template>
	<xsl:template match="*" mode="details_ontology">
		<xsl:variable name="ref">
			<xsl:choose>
				<xsl:when test="@rdf:ID">
					<xsl:value-of select="@rdf:ID" />
				</xsl:when>
				<xsl:when test="@rdf:about">
					<xsl:value-of select="@rdf:about" />
				</xsl:when>
				<xsl:otherwise>
					BLANK
				</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>

     <!--xsl:if test="string-length($ref)>0">	
		    <a name="{$ref}"><xsl:value-of select="$ref"/></a> 
     </xsl:if-->
		<table>
			<tr>
				<th colspan="2">
					<xsl:value-of select="//rdf:RDF/@xml:base" />
					<a>
						<xsl:value-of select="//owl:Ontology/rdfs:label" />
					</a>
					<span class="cp-type">
						(rdf:type
						<xsl:call-template name="url">
							<xsl:with-param name="ns" select="namespace-uri()" />
							<xsl:with-param name="name" select="local-name()" />
						</xsl:call-template>
						)
					</span>
				</th>
			</tr>
			<xsl:if test="count(*)+count(@*)>0">
				<xsl:apply-templates select="." mode="attribute" />
				<xsl:apply-templates select="*" mode="child" />
			</xsl:if>
		</table>
	</xsl:template>
	<xsl:template match="*" mode="details">
		<xsl:variable name="ref">
			<xsl:choose>
				<xsl:when test="@rdf:ID">
					<xsl:value-of select="@rdf:ID" />
				</xsl:when>
				<xsl:when test="@rdf:about">
					<xsl:value-of select="@rdf:about" />
				</xsl:when>
				<xsl:otherwise>
					BLANK
				</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>
		<table class="subsection" >
			<a name="{$ref}"/>
			<tbody>
				<xsl:if test="string-length($ref)>0">
					<tr>
						<th id="{substring-after(@rdf:about, '#')}">
							<xsl:choose>
								<xsl:when test="@rdf:ID">
									<a><xsl:value-of select="@rdf:ID" /></a>
								</xsl:when>
								<xsl:when test="@rdf:about">
									<xsl:call-template name="prettyUrl">
										<xsl:with-param name="name" select="@rdf:about"/>
									</xsl:call-template>
								</xsl:when>
								<xsl:otherwise>
									BLANK
								</xsl:otherwise>
							</xsl:choose>
							<span class="cp-type">
								(rdf:type
								<xsl:call-template name="url">
									<xsl:with-param name="ns" select="namespace-uri()" />
									<xsl:with-param name="name" select="local-name()" />
								</xsl:call-template>
								)
							</span>
						</th>
					</tr>
				</xsl:if>
				<xsl:if test="count(*)+count(@*)>0">
					<xsl:apply-templates select="." mode="attribute" />
					<xsl:apply-templates select="*" mode="child" />
				</xsl:if>
			</tbody>
		</table>
	</xsl:template>
	<xsl:template name="url">
		<xsl:param name="ns" />
		<xsl:param name="name" />
		<xsl:choose>
			<xsl:when test="$ns='&rdf;'">
				<a href="{concat($ns,$name)}">
					rdf:<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="$ns='&rdfs;'">
				<a href="{concat($ns,$name)}">
					rdfs:<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="$ns='&owl;'">
				<a href="{concat($ns,$name)}">
					owl:<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="$ns='&dc;'">
				<a href="{concat($ns,$name)}">
					dc:<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="$ns='&skos;'">
				<a href="{concat($ns,$name)}">
					skos:<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="$ns='&foaf;'">
				<a href="{concat($ns,$name)}">
					foaf:<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="$ns='&core;'">
				<a href="{concat($ns,$name)}">
					core:<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="$ns='&vg;'">
				<a href="{concat($ns,$name)}">
					vg:<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="$ns='http://purl.org/dc/terms/'">
                                <a href="{concat($ns,$name)}">
                                        dcterms:<xsl:value-of select="$name" />
                                </a>
                        </xsl:when>
			<xsl:when test="$ns='&spin;'">
                                <a href="{concat($ns,$name)}">
                                        spin:<xsl:value-of select="$name" />
                                </a>
                        </xsl:when>
			<xsl:when test="$ns=/rdf:RDF/@xml:base">
				<a href="{concat('#',$name)}">
					<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:when test="(string-length($ns)>0) or starts-with($name,'http://')">
				<a href="{concat($ns,$name)}">
					<xsl:value-of select="$name" />
				</a>
			</xsl:when>
			<xsl:otherwise>
				<xsl:value-of select="$name" />
			</xsl:otherwise>
		</xsl:choose>
	</xsl:template>
	<xsl:template name="prettyUrl">
		<xsl:param name="name" />
		<xsl:choose>
			<xsl:when test='starts-with($name, "&owl;")'>
				owl:<xsl:value-of select='substring-after($name, "&owl;")' />
			</xsl:when>
			<xsl:when test='starts-with($name, "&rdf;")'>
				rdf:<xsl:value-of select='substring-after($name, "&rdf;")' />
			</xsl:when>
			<xsl:when test='starts-with($name, "&rdfs;")'>
				rdfs:<xsl:value-of select='substring-after($name, "&rdfs;")' />
			</xsl:when>
			<xsl:when test='starts-with($name, "&xsd;")'>
				xsd:<xsl:value-of select='substring-after($name, "&xsd;")' />
			</xsl:when>
			<xsl:when test='starts-with($name, "&foaf;")'>
				foaf:<xsl:value-of select='substring-after($name, "&foaf;")' />
			</xsl:when>
			<xsl:when test='starts-with($name, "&skos;")'>
				skos:<xsl:value-of select='substring-after($name, "&skos;")' />
			</xsl:when>
			<xsl:when test='starts-with($name, "&core;")'>
				core:<xsl:value-of select='substring-after($name, "&core;")' />
			</xsl:when>
			<xsl:when test='starts-with($name, "&vg;")'>
				vg:<xsl:value-of select='substring-after($name, "&vg;")' />
			</xsl:when>
			<xsl:when test='starts-with($name,"http")'>
				<xsl:value-of select="$name" />
			</xsl:when>
			<xsl:when test="$name">
				<xsl:value-of select="$name" />
			</xsl:when>
		</xsl:choose>
	</xsl:template>
	<xsl:template match="*" mode="resource">
		<a href="{@rdf:resource}">
			<xsl:call-template name="prettyUrl">
				<xsl:with-param name="name" select="@rdf:resource"/>
			</xsl:call-template>
		</a>
	</xsl:template>
	<xsl:template match="*" mode="attribute">
		<xsl:for-each select="@*[local-name()!=lang]">
			<sup>
				<xsl:call-template name="url">
						<xsl:with-param name="ns" select="namespace-uri()" />
						<xsl:with-param name="name" select="local-name()" />
					</xsl:call-template>
					<xsl:call-template name="url">
						<xsl:with-param name="ns" select="''" />
						<xsl:with-param name="name" select="." />
					</xsl:call-template>
                        </sup>
		</xsl:for-each>
                 <xsl:if test="@rdf:datatype">
                 
		<xsl:for-each select="@rdf:datatype">
                   <sup>
                       <a href="{concat(namespace-uri(),.)}">
                        <xsl:call-template name="prettyUrl">
                          <xsl:with-param name="name" select="."/>
                        </xsl:call-template>
                      </a>
                     </sup>
		</xsl:for-each>
		</xsl:if>
	</xsl:template>
	<xsl:template match="*" mode="child">
		<tr>
			<td>
				<xsl:call-template name="url">
					<xsl:with-param name="ns" select="namespace-uri()" />
					<xsl:with-param name="name" select="local-name()" />
				</xsl:call-template>
			</td>
			<td>
				<xsl:choose>
					<xsl:when test="@rdf:resource">
						<xsl:apply-templates select="." mode="resource" />
					</xsl:when>
					<xsl:when test="local-name() = 'hasKey' and namespace-uri() = '&owl;'">
                                          <a href="{concat('&owl;', 'hasKey')}">owl:hasKey</a> consisting of (<br />
                                          <xsl:for-each select="descendant::*/rdf:first">
                                            <span class="indent">
                                              <xsl:apply-templates select="." mode="resource" />
                                            </span>
                                            <xsl:if test="position() != last()">
                                              <br /><em>and</em>
                                            </xsl:if>
                                            <br />
                                          </xsl:for-each>
                                          )
                                        </xsl:when>
					<xsl:when test="owl:Restriction">
						<xsl:apply-templates select="*" mode="restriction" />
					</xsl:when>
					<xsl:when test="@*">
						"<xsl:value-of select='text()' />"
						<xsl:if test="@xml:lang">@<xsl:apply-templates select="@xml:lang" mode="resource" />
						</xsl:if>
						<xsl:apply-templates select="." mode="attribute" />
					</xsl:when>
					<xsl:otherwise>
						<xsl:value-of select='text()' />
						<xsl:if test="@xml:lang">^^<xsl:apply-templates select="@xml:lang" mode="resource" />
						</xsl:if>
						<xsl:if test="owl:hasKey">
                                                  hasKey
                                                </xsl:if>
						<xsl:if test="owl:Restriction">
							<xsl:apply-templates select="*" mode="restriction" />
						</xsl:if>
						<xsl:if test="rdfs:Datatype" >
                                                    <xsl:apply-templates select="*" mode="restriction" />
                                                </xsl:if>
						<xsl:if test="owl:Class/owl:unionOf">
							<a href="http://www.w3.org/2002/07/owl#unionOf">owl:unionOf</a>
							<xsl:variable name="this">
								<xsl:value-of select="owl:Class/owl:unionOf/@rdf:nodeID" />
							</xsl:variable>
							(
							<xsl:choose>
								<xsl:when test="owl:Class/owl:unionOf[@rdf:parseType]/owl:Restriction">
									<xsl:for-each select="owl:Class/owl:unionOf[@rdf:parseType]">
                                                                                <xsl:if test="position() != last()">
                                                                                        <br />
                                                                                </xsl:if>
										<xsl:apply-templates select="*" mode="restriction" />
										<xsl:if test="position() != last()">
											<em>or</em><br />
										</xsl:if>
									</xsl:for-each>
								</xsl:when>
								<xsl:when test="owl:Class/owl:unionOf[@rdf:parseType]/owl:Class">
								<table>
									<xsl:for-each select="owl:Class/owl:unionOf[@rdf:parseType]">
										<xsl:for-each select='owl:Class'>
											<xsl:apply-templates select="*" mode="child" />
										</xsl:for-each>
									</xsl:for-each>
									</table>
								</xsl:when>
								<xsl:when test="owl:Class/owl:unionOf[@rdf:parseType]">
									<xsl:for-each select="owl:Class/owl:unionOf[@rdf:parseType]/child::*">
                                                                                <xsl:if test="position() = 1">
                                                                                        <br/>
                                                                                </xsl:if>
										<span class="indent"><a href="#{@rdf:about}">
											<xsl:value-of select="@rdf:about" />
										</a></span>
										<xsl:if test="position() != 1">
											<br /><em>or</em>
										</xsl:if>
										<br />
									</xsl:for-each>
								</xsl:when>
								<xsl:when test="owl:Class/owl:unionOf/rdf:Description">
                                                                        <xsl:for-each select="owl:Class/owl:unionOf/descendant::*/rdf:first">
                                                                                <xsl:if test="position() = 1">
                                                                                        <br/>
                                                                                </xsl:if>
                                                                                <span class="indent"><xsl:apply-templates select="." mode="resource" /></span>
                                                                                <xsl:if test="position() != last()">
                                                                                        <br /><em>or</em>
                                                                                </xsl:if>
                                                                                <br />
                                                                        </xsl:for-each>
                                                                </xsl:when>
								<xsl:otherwise>
									<xsl:for-each select="owl:Class/owl:unionOf/child::*">
                                                                                <xsl:if test="position() = 1">
                                                                                        <br/>
                                                                                </xsl:if>
										<span class="indent"><a href="{@rdf:about}">
											<xsl:value-of select="@rdf:about" />
										</a></span>
										<xsl:if test="position() != last()">
											<br /><em>or</em>
										</xsl:if>
										<br />
									</xsl:for-each>
								</xsl:otherwise>
							</xsl:choose>
							)
						</xsl:if>
						<xsl:if test="owl:Class/owl:intersectionOf">
							<a href="http://www.w3.org/2002/07/owl#intersectionOf">owl:intersectionOf</a>
							<xsl:variable name="this">
								<xsl:value-of select="owl:Class/owl:intersectionOf/@rdf:nodeID" />
							</xsl:variable>
							(
							<xsl:choose>
                                                                <xsl:when test="owl:Class/owl:intersectionOf/rdf:Description">
                                                                        <xsl:for-each select="owl:Class/owl:intersectionOf/descendant::*/rdf:first">
                                                                                <xsl:if test="position() != last()">
                                                                                        <br/>
                                                                                </xsl:if>
                                                                                <span class="indent"><xsl:apply-templates select="*" mode="restriction" /></span>
                                                                                <xsl:if test="position() != last()">
                                                                                        <em>and</em><br/>
                                                                                </xsl:if>
                                                                        </xsl:for-each>
                                                                </xsl:when>
								<xsl:when test="owl:Class/owl:intersectionOf[@rdf:parseType]">
									<xsl:for-each select="owl:Class/owl:intersectionOf[@rdf:parseType]">
                                                                                <xsl:if test="position() != last()">
                                                                                        <br/>
                                                                                </xsl:if>
										<span class="indent"><xsl:apply-templates select="*" mode="restriction" /></span>
										<xsl:if test="position() != last()">
                                                                                        <em>and</em><br/>
                                                                                </xsl:if>
									</xsl:for-each>
								</xsl:when>
								<xsl:when test="owl:Class/owl:intersectionOf">
									<xsl:for-each select="owl:Class/owl:unionOf/child::*">
                                                                                <xsl:if test="position() != last()">
                                                                                        <br/>
                                                                                </xsl:if>
										<span class="indent"><a href="#{@rdf:about}">
											<xsl:value-of select="@rdf:about" />
										</a></span>
										<xsl:if test="position() != last()">
                                                                                        <em>and</em><br/>
                                                                                </xsl:if>
									</xsl:for-each>
								</xsl:when>
							</xsl:choose>
							)
						</xsl:if>
					</xsl:otherwise>
				</xsl:choose>
			</td>
		</tr>
	</xsl:template>
	<xsl:template match="*" mode="restriction">
		Restrict
		<a href="#{owl:onProperty/@rdf:resource}">
                      <xsl:call-template name="prettyUrl">
                         <xsl:with-param name="name" select="owl:onProperty/@rdf:resource"/>
                      </xsl:call-template>
		</a> 
		<xsl:if test="owl:someValuesFrom/@rdf:resource">
			to <a href="{concat('&owl;', 'someValuesFrom')}">owl:someValuesFrom </a>
			<xsl:apply-templates select="owl:someValuesFrom" mode="resource"/>
		</xsl:if>
		<xsl:if test="owl:allValuesFrom/@rdf:resource">
			to <a href="{concat('&owl;', 'allValuesFrom')}">owl:allValuesFrom </a>
			<xsl:apply-templates select="owl:allValuesFrom" mode="resource"/>
		</xsl:if>
		<xsl:if test="owl:cardinality">
			to <a href="{concat('&owl;', 'cardinality')}">owl:cardinality </a>
				<a href="#{owl:cardinality/@rdf:datatype}">
				<xsl:value-of select="owl:cardinality/text()" />
			</a>
		</xsl:if>
		<xsl:if test="owl:maxCardinality">
			to <a href="{concat('&owl;', 'maxCardinality')}">owl:maxCardinality </a>
			<a href="#{owl:maxCardinality/@rdf:datatype}">
				<xsl:value-of select="owl:maxCardinality/text()" />
			</a>
		</xsl:if>
		<xsl:if test="owl:minCardinality">
			to <a href="{concat('&owl;', 'minCardinality')}">owl:minCardinality </a>
			<a href="#{owl:minCardinality/@rdf:datatype}">
				<xsl:value-of select="owl:minCardinality/text()" />
			</a>
		</xsl:if>
		<xsl:if test="owl:onClass">
			<a href="{concat('&owl;', 'onClass')}"> owl:onClass </a>
			<xsl:apply-templates select="owl:onClass" mode="resource"/>
		</xsl:if>
		<xsl:if test="owl:onDataRange">
			<a href="{concat('&owl;', 'onDataRange')}"> owl:onDataRange </a>
			<xsl:apply-templates select="owl:onDataRange" mode="resource"/>
		</xsl:if>
		<xsl:if test="owl:onDatatype">
                        <a href="{concat('&owl;', 'onDataType')}"> owl:onDatatype </a> to 
                        <xsl:apply-templates select="owl:onDatatype" mode="resource"/> with a value of 
                        <xsl:for-each select="owl:onDatatype/../owl:withRestrictions/rdf:Description/descendant::*/rdf:Description/child::*">
                          <xsl:if test="local-name() = 'minInclusive'" >
                            at least <a href="#{@rdf:datatype}" ><xsl:value-of select="text()" /> </a>
                          </xsl:if>
                          <xsl:if test="local-name() = 'maxInclusive'" >
                            at most <a href="#{@rdf:datatype}" ><xsl:value-of select="text()" /> </a>
                          </xsl:if>
                        </xsl:for-each>
                </xsl:if>
		<xsl:if test="owl:maxQualifiedCardinality">
			to <a href="{concat('&owl;', 'maxQualifiedCardinality')}">owl:maxQualifiedCardinality </a>
			<a href="#{owl:maxQualifiedCardinality/@rdf:datatype}">
                            <xsl:value-of select="owl:maxQualifiedCardinality/text()" />
			</a>
		</xsl:if>
		<xsl:if test="owl:minQualifiedCardinality">
			to <a href="{concat('&owl;', 'minQualifiedCardinality')}">owl:minQualifiedCardinality </a>
			<a href="#{owl:minQualifiedCardinality/@rdf:datatype}">
				<xsl:value-of select="owl:minQualifiedCardinality/text()" />
			</a>
		</xsl:if>
		<xsl:if test="owl:hasValue">
			<a href="{concat('&owl;', 'hasValue')}"> owl:HasValue </a>
			<a href="#{owl:hasValue/@rdf:datatype}">
				<xsl:value-of select="owl:hasValue/text()" />
			</a>
		</xsl:if>
		<br/>
	</xsl:template>
</xsl:stylesheet>
