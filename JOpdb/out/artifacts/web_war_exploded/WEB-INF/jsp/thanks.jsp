<%--
  Created by IntelliJ IDEA.
  User: lyd
  Date: 10/17/17
  Time: 11:39 AM
  To change this template use File | Settings | File Templates.
--%>
<%--
  Created by IntelliJ IDEA.
  User: lyd
  Date: 10/13/17
  Time: 7:12 PM
  To change this template use File | Settings | File Templates.
--%>
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
    <title>OPDB | thanks</title>
</head>
<body class="loading">
<div id="wrapper">
    <div id="bg"></div>
    <div id="overlay"></div>
    <div id="main">

        <!-- Header -->
        <header id="header">
            <h1>OPDB</h1>
            <nav><hr/>
                <p>Thanks to help us complete our database!</p>
                <p>Please check your email inbox!</p>
            </nav>
        </header>

        <!-- Footer -->
        <footer id="footer">
            <span class="copyright">&copy; github: <a href="https://github.com/GodInLove/OPDB">OPDB</a>.</span>
        </footer>

    </div>
</div>
<!--[if lte IE 8]><script src="assets/js/ie/respond.min.js"></script><![endif]-->
<script>
    window.onload = function() { document.body.className = ''; }
    window.ontouchmove = function() { return false; }
    window.onorientationchange = function() { document.body.scrollTop = 0; }
</script>
</body>
</html>