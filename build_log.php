<html>
<head>
  <link href="http://www.dealii.org/screen.css" type="text/css" rel="StyleSheet">
  <title>deal.II Build Test Results</title>
</head>
<body>
<h1>Build Test Results</h1>

<?php
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

  echo "<table border=1>\n<tr>";
  $n = pg_num_fields($result);
  for ($i=0;$i<$n;++$i)
  {
    echo $th, pg_field_name($result, $i), "</th>";
  }
  echo "</tr>\n";

  while ($row = pg_fetch_row($result))
  {
    echo "<tr>";
    for ($i=0; $i<$n; ++$i)
    {
      echo "$td $row[$i] </td>";
    }
    echo "</tr>\n";
  }
  echo "</table>\n";
}


function view_log ($connectstring, $id)
{
  $database = pg_pconnect("$connectstring");
  
  if (!$database)
  {
    echo "Connection to database failed.";
    exit;
  }
  $query = "select log from log where id=$id;";
  $result = pg_query($database, $query);
  if (!$result)
  {
    echo "Could not execute query: $query";
    exit;
  }
  
  $row = pg_fetch_row($result);
  echo "<pre>\n$row[0]\n</pre>\n";
}
?>


<?php
  import_request_variables("pg", "form_");
  echo "<h2>Details for configuration $form_id</h2>";
  view_table("dbname=deal", "build_history", "where id=$form_id", "build");
  view_log("dbname=deal", $form_id);
?>
</table>

<h2>Previous Test Page</h2>

The old page with build test results can be found <a
       href="http://www.dealii.org/cgi-bin/show_build.pl?summary=1">here</a>.

</body>
</html>
