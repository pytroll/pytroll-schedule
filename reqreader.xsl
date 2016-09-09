<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
	<xsl:output method="html" indent="yes" encoding="UTF-8" />

	<xsl:template match="acquisition-schedule">
		<html>
			<head>
				<meta charset="UTF-8" />
				<title>Überflugplan
				</title>
				<script type="text/javascript">
					function switcher(src) {
					bild = document.getElementById("plot");
					bild.src = src;
					bild.visible = true;
					}
				</script>
			</head>
			<body>
				<h1>Überflugplan</h1>
				<h2>
					Erstellt von:
					<xsl:value-of select="properties/requested-by" />
					um:
					<xsl:value-of select="properties/requested-on" />
				</h2>
				<h2>
					Antenne:
					<xsl:value-of select="properties/station" />
				</h2>
				<table width="90%">
					<tr>
						<td>
							<table style="border-spacing:12px;">
								<tr>
									<th>Pass</th>
									<th>Satellit</th>
									<th>Aufgang</th>
									<th>Untergang</th>
								</tr>
								<xsl:apply-templates select="child::node()" />
							</table>
						</td>
						<td valign="top" align="right" width="30%">
							<img id="plot" src="" width="400px" height="300px" visible="false" />
						</td>
					</tr>
				</table>
			</body>
		</html>
	</xsl:template>

	<xsl:attribute-set name="css">
		<xsl:attribute name="style">font-size:10pt;</xsl:attribute>
	</xsl:attribute-set>

	<xsl:template match="pass">
		<xsl:variable name="plotfullname">
			<xsl:value-of select="./@img" />
		</xsl:variable>
		<xsl:variable name="plotname">
			<xsl:value-of select="substring-after($plotfullname, 'plots.ofb' )" />
		</xsl:variable>
		<xsl:element name="tr">
			<xsl:attribute name="onClick">javascript:switcher('plots.ofb<xsl:value-of select="$plotname" />')</xsl:attribute>
			<xsl:choose>
				<xsl:when test="./@rec = 'True'">
					<xsl:attribute name="style">font-family:monospace; font-weight:bold;</xsl:attribute>
				</xsl:when>
				<xsl:otherwise>
					<xsl:attribute name="style">font-family:monospace; font-weight:normal;</xsl:attribute>
				</xsl:otherwise>
			</xsl:choose>
			<xsl:call-template name="passvalues">
				<xsl:with-param name="sat" select="." />
			</xsl:call-template>
		</xsl:element>
	</xsl:template>

	<xsl:template name="passvalues">
		<xsl:param name="sat" />
		<td>
			<xsl:number format="1" />
		</td>
		<xsl:element name="td" use-attribute-sets="css">
			<xsl:value-of select="./@satellite" />
		</xsl:element>
		<xsl:element name="td" use-attribute-sets="css">
			<xsl:value-of select="./@start-time" />
		</xsl:element>
		<xsl:element name="td" use-attribute-sets="css">
			<xsl:value-of select="./@end-time" />
		</xsl:element>
	</xsl:template>

	<xsl:template match="*" />
</xsl:stylesheet>
