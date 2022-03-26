<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
 
<xsl:template match="/">
<html>
<body>
  <h2>Logs</h2>
  <table frame="none" cellpadding="3">
    <xsl:for-each select="log/logentry">
      <tr style="background:#FFFFFF">
	<td><xsl:value-of select="@revision"/></td>
	<td><xsl:value-of select="substring(date, 1, 10)"/></td>
	<td><xsl:value-of select="author"/></td>
	<td><xsl:value-of select="msg"/></td>
      </tr>
      <xsl:for-each select="paths/path">
	<tr style="background:#FFFFFF">
	  <td></td>
	  <td></td>
	  <td style="background:#AAAAAA"><xsl:value-of select="@action"/></td>
	  <td style="background:#AAAAAA"><xsl:value-of select="."/></td>
	</tr>
      </xsl:for-each>
    </xsl:for-each>
  </table>
</body>
</html>
</xsl:template>

</xsl:stylesheet>
