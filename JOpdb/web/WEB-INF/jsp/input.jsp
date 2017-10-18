<%--
  Created by IntelliJ IDEA.
  User: zhaoliang
  Date: 2017/10/13
  Time: 下午5:53
  To change this template use File | Settings | File Templates.
--%>
<!DOCTYPE HTML>
<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<html>
<head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <!--[if lte IE 8]><script src="../../static/js/ie/html5shiv.js"></script><![endif]-->
    <link rel="stylesheet" href="../../static/css/main.css" />
    <!--[if lte IE 8]><link rel="stylesheet" href="../../static/css/ie8.css" /><![endif]-->
    <!--[if lte IE 9]><link rel="stylesheet" href="../../static/css/ie9.css" /><![endif]-->
    <title>OPDB | predict</title>
</head>
<body class="loading">
<div id="wrapper">
    <div id="bg"></div>
    <div id="overlay"></div>
    <div id="main">

        <!-- Header -->
        <header id="header">
            <h1>OPDB</h1>
            <p>SCU &nbsp;&bull;&nbsp; lyd &nbsp;&bull;&nbsp; operon prediction database</p>
            <nav>
                <form action="/operon/post.do" method="post" class="">
                    <input id="email" type="text" placeholder="your@email.com" required id = "email" name="emailtext"><br>
                    <hr/>
                    <input id="kegg_id" type="text" placeholder="eco" required id ="kegg_id" name="kegg_idtext"><br>
                    <input id="srr_num" type="text" placeholder="SRR5486953" required id ="srr_num" name="srr_numtext"><br>
                    <select id="method" required id = "method" name="methodtext">
                        <option value="rockhopper">rockhopper</option>
                    </select><br>
                    <input type="submit" value="submit"/>

                </form>
            </nav>
            <p>We will send you a email if the job finished.</p>
        </header>

        <!-- Footer -->
        <footer id="footer">
            <span class="copyright">&copy; github: <a href="https://github.com/GodInLove/OPDB">OPDB</a>.</span>
        </footer>

    </div>
</div>
<!--[if lte IE 8]><script src="../../static/js/ie/respond.min.js"></script><![endif]-->
<script>
    window.onload = function() { document.body.className = ''; }
    window.ontouchmove = function() { return false; }
    window.onorientationchange = function() { document.body.scrollTop = 0; }
</script>
</body>
</html>