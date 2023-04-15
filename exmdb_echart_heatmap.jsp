<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
        <meta charset="utf-8">
    </head>
    <body style="height: 100%; margin: 0">
        <div id="container" style="height: 800px"></div>

        
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5/dist/echarts.min.js"></script>
        <!-- Uncomment this line if you want to dataTool extension
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5/dist/extension/dataTool.min.js"></script>
        -->
        <!-- Uncomment this line if you want to use gl extension
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts-gl@2/dist/echarts-gl.min.js"></script>
        -->
        <!-- Uncomment this line if you want to echarts-stat extension
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts-stat@latest/dist/ecStat.min.js"></script>
        -->
        <!-- Uncomment this line if you want to use map
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5/map/js/china.js"></script>
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5/map/js/world.js"></script>
        -->
        <!-- Uncomment these two lines if you want to use bmap extension
        <script type="text/javascript" src="https://api.map.baidu.com/api?v=2.0&ak=<Your Key Here>"></script>
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@{{version}}/dist/extension/bmap.min.js"></script>
        -->

        <script type="text/javascript">
var dom = document.getElementById("container");
var myChart = echarts.init(dom);
var app = {};

var option;

<%
String cancername = "GSE142213";  
String path = request.getRealPath("/").replace("\\","/");   //get root real path 
String heatmap  = path+"echart_heatmap/"+cancername+"_heatmap.txt";

BufferedReader br = new fcon(heatmap).getBr();
String str = br.readLine(); 

while(str!=null && str.length()>0) {
	out.println(str+";\n");
	str = br.readLine();
}

%>

// prettier-ignore

option = {
  tooltip: {
    position: 'top'
  },
  grid: {
    height: '100%',
    top: '10%',
    containLabel: true
  },
  xAxis: {
      type: 'category',
      data: rowname,
      splitArea: {
          show: true
      },
     axisLabel: {  
     interval:0,  
     rotate:90,
     textStyle: {
              color: '#000',
              fontSize:'5',
              itemSize:''
     }
  }  
  },
  yAxis: {
      type: 'category',
      data: colname,
      splitArea: {
          show: true
      },
      axisLabel: {  
      interval:0,  
      textStyle: {
              color: '#000',
              fontSize:'5',
              itemSize:''
      }
  }
  },
  visualMap: {
    min: 0,
    max: 3,
    calculable: true,
    orient: 'horizontal',
    left: 'center',
    bottom: '0%'
  },
  series: [
    {
      name: 'Punch Card',
      type: 'heatmap',
      data: data,
      label: {
        show: false
      },
      emphasis: {
        itemStyle: {
          shadowBlur: 10,
          shadowColor: 'rgba(0, 0, 0, 0.5)'
        }
      }
    }
  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

        </script>
    </body>
</html>
    