<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <!-- Formatting aquisition-schedule-report_ofb.xml for nice display in browser. -->
  <!-- DWD/amaul 13.09.2016 -->
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
          bild.style.visibility = "visible";
          }
        </script>
      </head>
      <body>
        <div style="position:fixed; height:15%; ">
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
        </div>
        <div style="position:absolute; height:85%; bottom:0px; overflow:auto; width:45%; margin-left:2%; ">
          <table style="border-spacing:12px;">
            <tr>
              <th>Pass</th>
              <th>Satellit</th>
              <th>Aufnahme-Anfang</th>
              <th>Aufnahme-Ende</th>
            </tr>
            <xsl:apply-templates select="child::node()" />
          </table>
        </div>
        <div style="position:fixed; width:50%; top:300px; right:0px; text-align:center;">
          (Click auf Überflug zeigt Plot)
        </div>
        <div style="position:fixed; width:50%; top:200px; right:0px; text-align:center; ">
          <img id="plot" src="" style="width:600px; height:450px; visibility:hidden;" />
        </div>
      </body>
    </html>
  </xsl:template>

  <xsl:attribute-set name="css">
    <xsl:attribute name="style">font-size:10pt;</xsl:attribute>
  </xsl:attribute-set>

  <xsl:template match="pass">
    <xsl:variable name="plotname">
      <xsl:value-of select="substring-after(./@img, 'plots.ofb' )" />
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
