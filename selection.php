<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
       "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
  <title>deal.II Homepage Gallery Selection</title>
  <link href="screen.css" type="text/css" rel="StyleSheet">
  <meta name="author" content="the deal.II authors <authors@dealii.org>">
  <meta name="copyright" content="Copyright (C) 2003, 2004, 2005, 2007, 2010, 2012, 2013 by the deal.II authors">
  <meta name="date" content="$Date: 2013-08-22 02:51:54 +0200 (Thu, 22 Aug 2013) $">
  <meta name="svn_id" content="$Id: selection.php 4952 2013-08-22 00:51:54Z ted.studley $">
  <meta name="robots" content="noindex">

  <script type="text/javascript" src="js/jquery.min.js"></script>
  <script type="text/javascript" src="js/tinycarousel.js"></script>
  <script type="text/javascript">
    $(document).ready(function() {
      $('#slider-code').tinycarousel({});
    });
  </script>
</head>

<body class="navbar">

<p align="center"><img src="pictures/logo100.png" alt="logo"></p>
<hr>

  <b><small>Selected images</small></b>

<?php
  srand((float) microtime() * 10000000);
  $d = dir("images/wiki/gallery/thumbs");
  $i = 0;
  $list = array();
  while (false !== ($entry = $d->read()))
  {
    if (preg_match("/\.[a-z]+/", $entry))
    {
      $text = "<p align=\"center\"><a target=\"_top\""
	." href=\"images/wiki/gallery/"
        .$entry."\"><img width=\"100px\" src=\"images/wiki/gallery/thumbs/"
        .$entry."\" border=\"0\" alt=\"".$entry."\"></a></p>\n";
      array_push($list, $text);
    } else {
      print "<!-- $entry -->\n";
    }
  }
  $d->close();
?>

<?php
  $out = array();
  $out = array_rand($list, 8);
  foreach ($out as $entry)
  {
    print $list[$entry];
  }
?>

<hr>

</body>
</html>
