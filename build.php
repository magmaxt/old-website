<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
                 "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
  <link href="http://www.dealii.org/screen.css" type="text/css" rel="StyleSheet">
  <title>deal.II Build Test Results</title>
</head>
<body>

<?php
function overview ($connectstring, $table, $conditions, $css_class)
{
  $database = pg_pconnect("$connectstring");
  
  if (!$database)
  {
    echo "Connection to database failed.";
    exit;
  }
  $query = "select * from $table $conditions;";
  $result = pg_query($database, $query);
  if (!$result)
  {
    echo "Could not execute query: $query";
    exit;
  }  

  $n = pg_num_fields($result) - 1;
  $n_days = 15;

  // build an array of past days that holds a date id and a string representing it
  // we will use the id to later sort test results into the correct column
  for ($i=0; $i<$n_days; ++$i)
    {
      $past_date[$i] = mktime(0, 0, 0, date("m")  , date("d")-$i, date("Y"));
      $past_date_string[$i] = date("Y/<\\b\\r>m/<\\b\\r>d", $past_date[$i]);
    }

  // write header of table
  echo "  <tr>\n";
  echo "    <td class=\"$css_class\">Date <br> (CE(S)T timezone)</td>\n";
  for ($i=0; $i<$n_days; ++$i)
    {
      echo "    <td class=\"$css_class\">$past_date_string[$i]</td>\n";
    }
  echo "  </tr>\n";

  // now loop over all configurations
  while ($row = pg_fetch_row($result))
  {
    // start with an empty row of test results
    for ($i=0; $i<$n_days; ++$i)
      {
         $past_results[$i] = "";
      }

    // then read all test results for this configuration and loop over them
    $history_query = "select * from build_history where id=$row[$n];";
    $history_result = pg_query($database, $history_query);

    while ($history_row = pg_fetch_row($history_result))
      {
        // determine the date of this build result, and add it to the
	// previous results from this date
        list($date,$time) = split(" ", $history_row[1]);
	list($year,$month,$day) = split("-", $date);
	$this_date = mktime(0, 0, 0, $month, $day, $year);
        for ($i=0; $i<$n_days; ++$i)
          {
            if ($this_date == $past_date[$i])
              {
                $past_results[$i] = "$past_results[$i]$history_row[0]";
              }
          }
      }

    // check if this configuration has produced anything within the last $n_days
    // days; if not, don't produce output
    $has_something = 0;
    for ($i=0; $i<$n_days; ++$i)
      {
        if ($past_results[$i] != "")
          {
            $has_something = 1;
            break;
          }
      }
    if ($has_something == 0)
      {
        continue;
      }

    // if this configuration did have results recently, output a row of data
    echo "  <tr>\n";
    echo "    <td>\n";
    echo "      <a href=\"#configuration-$row[$n]\">History for $row[$n]</a>\n";
    echo "    </td>\n";

    for ($i=0; $i<$n_days; ++$i)
      {
        echo "      <td class=\"$css_class\">\n";
        echo "        $past_results[$i]\n";
	echo "      </td>\n";
      }
    echo "  </tr>\n";
  }
}


function view_table ($connectstring, $table, $conditions, $css_class)
{
  $database = pg_pconnect("$connectstring");
  if ($css_class != "")
  {
    $th = "<th class=\"$css_class\">";
    $td = "<td class=\"$css_class\">";
  } else {
    $th = "<th>";
    $td = "<td>";
  }
  
  if (!$database)
  {
    echo "Connection to database failed.";
    exit;
  }
  $query = "select * from $table $conditions;";
  $result = pg_query($database, $query);
  if (!$result)
  {
    echo "Could not execute query: $query";
    exit;
  }  

  $n = pg_num_fields($result) - 1;

  while ($row = pg_fetch_row($result))
  {
    echo "  <tr>\n";

    for ($i=0; $i<$n; ++$i)
    {
      $text = $row[$i];
      $text = ereg_replace("x86_64-unknown-linux-gnu", "Linux64", $text);
      $text = ereg_replace("i686-pc-linux-gnu", "Linux", $text);
      $text = ereg_replace("sparc-sun-s", "S", $text);
      $text = ereg_replace("i686-pc-cygwin", "Windows Cygwin", $text);
      $text = ereg_replace("powerpc-apple-darwin", "Apple Darwin", $text);

      $text = ereg_replace("Rannacher-Sekretariat-Computer-OSX", "McChef", $text);

      echo "    $td $text </td>\n";
    }
    echo "    <td>\n";
    echo "      <a name=\"configuration-$row[$n]\"></a>\n";
    echo "      <form action=\"build_log.php\" method=\"post\">\n";
    echo "        <input type=\"hidden\" name=\"id\" value=\"$row[$n]\">\n";
    echo "        <input type=\"submit\" value=\"Details for $row[$n]\">\n";
    echo "      </form>\n";
    echo "    </td>\n";
    echo "  </tr>\n";
  }
}
?> 


<h1>Build Test Results</h1>

<h2>Overview of recent build test results</h2>

<table border=1 align="center">
<?php 
  overview("dbname=deal", "build_now", "order by id desc", "build");
?>
</table>


<h2>Details of configurations</h2>


<table border=1>
  <tr>
    <th></th>
    <th colspan="3">Configuration</th>
    <th colspan="3">Options</th>
    <th colspan="5">Library support</th>
  </tr>

  <tr>
    <th>OK</th>
    <th><a href="?order=compiler">compiler</a></th>
    <th><a href="?order=system">system</a></th>
    <th>version</th>
    <th>shared</th>
    <th>fparser</th>
    <th>MT</th>
    <th>blas</th>
    <th>lapack</th>
    <th>metis</th>
    <th>petsc</th>
    <th>umfpack</th>
    <th><a href="build.php">age</a></th>
    <th><a href="build.php?order=host">Host</a></th>
    <th><a href="build.php?order=id">Id</a></th>
  </tr>

<?php
    if ($_GET['order'] == 'compiler')
    {
      view_table("dbname=deal", "build_now", "order by compiler DESC, system, version DESC", "build");
    } else if ($_GET['order'] == 'system') {
      view_table("dbname=deal", "build_now", "order by system, compiler DESC, version DESC", "build");
    } else if ($_GET['order'] == 'id') {
      view_table("dbname=deal", "build_now", "order by id DESC", "build");
    } else if ($_GET['order'] == 'host') {
      view_table("dbname=deal", "build_now", "order by host DESC", "build");
    } else {
      view_table("dbname=deal", "build_now", "", "build");
    }
  ?>
</table>
<div class="right">
  <p>
    <a href="http://validator.w3.org/check?uri=referer" target="_top"><img
        style="border:0"
        src="http://www.w3.org/Icons/valid-html401"
	alt="Valid HTML 4.01!"></a>
    <a href="http://jigsaw.w3.org/css-validator/check/referer" target="_top"><img
       style="border:0;width:88px;height:31px"
       src="http://jigsaw.w3.org/css-validator/images/vcss" 
       alt="Valid CSS!">
 </a>
 </p>
</div>
</body>
</html>
