<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx" style="height: 100%">
    <head>
        <meta charset="utf-8">
    </head>
    <body style="height: 100%; margin: 0">
        <div id="container" style="height: 100%"></div>

        
        <script type="text/javascript" src="echart/echarts.min.js"></script>



        <script type="text/javascript">
			var dom = document.getElementById("container");
			var myChart = echarts.init(dom);
			var app = {};
			
			var option;
			
			<% 
			
			String searchname = "AC007728.2";  
			if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){
				searchname = request.getParameter("searchname");
			}
			

			String gseid = "COVID19_GSE154567_recovery_count";  
			if(request.getParameter("gseid")!=null&&!request.getParameter("gseid").equals("")){
				gseid = request.getParameter("gseid");
			}

			
			double QHT_array[] = new double[64];
			String exp = dbhello.exmdb_get_per_exp(searchname,gseid);
			
			
			
			if(exp.equals("")){
// 				System.out.println("NULL exp: "+ exp);
				exp = "#";
				searchname = "No expression in this data";
			}

			String gene_index [] = exp.split("#");
				for(int i=0;i<gene_index.length;i++){
					String index_exp[] =  gene_index[i].split("_");
//	 				System.out.println(Integer.parseInt(index_exp[0]));
// 	 				System.out.println(gene_index[i]);
					QHT_array[(Integer.parseInt(index_exp[0])-1)] = Double.parseDouble(index_exp[1]);
				}

			
			
// 			System.out.println(Arrays.toString(QHT_array));
			String qht = dbhello.getpercent(QHT_array);
// 			System.out.println(qht);
			out.println(qht);
			%>
			
			option = {
			              title: {
			                  text: '<%=searchname%>',
			                  textStyle:{fontFamily: 'Montserrat'},
			              },
			  tooltip: {
			    trigger: 'axis',
			    textStyle:{fontFamily: 'Montserrat',fontSize:15},
			    axisPointer: {
			      // Use axis to trigger tooltip
			      type: 'shadow' // 'shadow' as default; can also be 'line' or 'shadow'
			    }
			  },
			  color:["#F4B71D","#114182"],
			  legend: {
				  textStyle:{fontFamily: 'Montserrat',fontSize:15},
				  left: 'center',
				  type: 'scroll'
			  },
			  grid: {
			    left: '3%',
			    right: '4%',
			    bottom: '3%',
			    containLabel: true
			  },
			  xAxis: {
			      axisLabel: {
			          fontFamily: 'Montserrat',
			          rotate: 45,
			      },
			    type: 'category',
			    data: label
			    
			  },
			  yAxis: {
			      axisLabel: {
			          rotate: 45,
			          fontFamily: 'Montserrat'
			      },
			    type: 'value'
			  },
			  series: [
			
			    {
			      name: 'Exp',
			      type: 'bar',
			      stack: 'total',
			      label: {
			        show: true,
			        fontFamily: 'Montserrat'
			      },
			      emphasis: {
			        focus: 'series'
			      },
			      data: exp
			    },
			    {
			      name: 'No Exp',
			      type: 'bar',
			      stack: 'total',
			      label: {
			        show: true,
			        fontFamily: 'Montserrat'
			      },
			      emphasis: {
			        focus: 'series'
			      },
			      data: Noexp
			    },
			  ]
			};
			
			if (option && typeof option === 'object') {
			    myChart.setOption(option);
			}

        </script>
    </body>
</html>
    