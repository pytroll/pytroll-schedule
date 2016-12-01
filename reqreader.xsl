<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <!-- Formatting aquisition-schedule-report_ofb.xml for nice display in browser. -->
  <!-- DWD/amaul 13.09.2016 -->
  <xsl:output method="html" indent="yes" encoding="UTF-8" />

  <xsl:template match="acquisition-schedule">
    <html>
      <head>
        <meta charset="UTF-8" />
        <title>Überflugplan</title>
        <style>
          body { font-family:sans-serif; font-size:10pt; font-weight:normal; }
          h1 { font-size:16pt; font-weight:bold; }
          h2 { font-size:12pt; font-weight:bold; }
          div.header { position:fixed; max-height:15%; }
          div.passes { position:absolute; height:85%; bottom:0px; overflow:auto; width:45%; margin-left:2%; font-family:monospace; }
          div.plot { position:fixed; width:50%; top:200px; right:0px; text-align:center; }
          div.plot_rmk { position:fixed; width:50%; top:300px; right:0px; text-align:center; }
          img.plot_img { width:600px; height:450px; visibility:hidden; }
          .rec_true { font-weight:bold; }
          .rec_false { font-weight:normal; }
        </style>
        <script type="text/javascript">
          function switcher(src) {
          bild = document.getElementById("plot");
          bild.src = src;
          bild.style.visibility = "visible";
          }
        </script>
      </head>
      <body>
        <div class="header">
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
        <div class="passes">
          <table style="border-spacing:12px;">
            <tr>
              <th>Pass</th>
              <th>Satellit</th>
              <th>Datum</th>
              <th>Anfang</th>
              <th>Ende</th>
            </tr>
            <xsl:apply-templates select="child::node()" />
          </table>
        </div>
        <div class="plot_rmk">
          (Click auf Überflug zeigt Plot)
        </div>
        <div class="plot">
          <img id="plot" src="" class="plot_img" />
        </div>
      </body>
    </html>
  </xsl:template>

  <xsl:template match="pass">
    <xsl:variable name="plotname">
      <xsl:value-of select="substring-after(./@img, 'plots.' )" />
    </xsl:variable>
    <xsl:element name="tr">
      <xsl:attribute name="onClick">javascript:switcher('plots.<xsl:value-of select="$plotname" />')</xsl:attribute>
      <xsl:choose>
        <xsl:when test="./@rec = 'True'">
          <xsl:attribute name="class">rec_true</xsl:attribute>
        </xsl:when>
        <xsl:otherwise>
          <xsl:attribute name="class">rec_false</xsl:attribute>
        </xsl:otherwise>
      </xsl:choose>
      <xsl:call-template name="passvalues">
        <xsl:with-param name="sat" select="." />
      </xsl:call-template>
    </xsl:element>
  </xsl:template>

  <xsl:template name="passvalues">
    <xsl:param name="sat" />
    <xsl:element name="td">
      <xsl:number format="1" />
    </xsl:element>
    <xsl:element name="td">
      <xsl:value-of select="./@satellite" />
    </xsl:element>
    <xsl:element name="td">
      <xsl:value-of select="substring(@start-time, 0, 11)" />
    </xsl:element>
    <xsl:element name="td">
      <xsl:value-of select="substring(@start-time, 12)" />
    </xsl:element>
    <xsl:element name="td">
      <xsl:value-of select="substring(@end-time, 12)" />
    </xsl:element>
  </xsl:template>

  <xsl:template match="*" />
</xsl:stylesheet>
