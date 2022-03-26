<html>
<head>
  <title>deal.II Homepage Gallery Selection</title>
  <link href="screen.css" type="text/css" rel="stylesheet">
  <link href="css/sidebar-carousel2.css" type="text/css" rel="stylesheet">
  <meta name="author" content="the deal.II authors <authors@dealii.org>">
  <meta name="copyright" content="Copyright (C) 2003, 2004, 2005, 2007, 2010, 2012, 2013, 2014 by the deal.II authors">
  <meta name="date" content="$Date:$">
  <meta name="svn_id" content-"$Id:$">
  <meta name="robots" content="noindex">

  <script type="text/javascript" src="js/jquery.min.js"></script>
  <script type="text/javascript" src="js/jquery.tinycarousel.min.js"></script>
  <script type="text/javascript">
    $(document).ready(function() {
      $('#slider-code').tinycarousel({ axis: 'y', interval: true, display: 1, intervaltime: 4000 });
    });
  </script>
</head>

<body class="navbar">

  <p align="center"><img src="pictures/logo100.png" alt="logl"></p>
  <hr>
  <?php
    $d = dir("images/wiki/gallery");
    $i = 0;
    $list = array();
    while (false !== ($entry = $d->read()))
    {
      if (preg_match("/\.[a-z]+/", $entry))
      {
        $text = "<li><a target=\"_top\""
          ." href=\"images/wiki/gallery/"
          .$entry."\"><img width='100%' src=\"images/wiki/gallery/"
          .$entry."\" alt=\"".$entry."\"></a></li>\n";
        array_push($list, $text);
      }
    }
    $d->close();
  ?>
  <div id="slider-code">
    <div class="viewport">
      <ul class="overview">
        <?php
          foreach ($list as $value) 
          {
            print $value;
          }
        ?>
      </ul>
    </div>
  </div>
