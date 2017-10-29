package edu.scu.controller;

import edu.scu.entity.Operon;
import edu.scu.service.MailService;
import edu.scu.service.OperonService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpSession;
import java.io.File;
import java.io.IOException;
import java.util.List;

@Controller
@RequestMapping("/operon")
public class OperonController {
    private final OperonService operonService;
    private final MailService mailService;

    @Autowired
    public OperonController(OperonService operonService, MailService mailService) {
        this.operonService = operonService;
        this.mailService = mailService;
    }

    @RequestMapping(value = "/post.do")
    public String postParameter(HttpServletRequest request,Model model) {
        String sentTo = request.getParameter("emailtext");
        String kegg_id = request.getParameter("kegg_idtext");
        String srr_num = request.getParameter("srr_numtext");
        String method = request.getParameter("methodtext");

        String isNeed = operonService.isNeed(srr_num);
        if (isNeed.equals("need")){
            String predictResult = operonService.predictOP(srr_num, kegg_id , method);
            if (predictResult.equals("done")){
                mailService.sendMail(sentTo,
                        "OPDB job",
                        "your job link is http://bioinfor.scu.edu.cn/OPDB/JBrowse/index.html?data=" + srr_num + "\n\ncontact us with ttttttliu@qq.com");
            }
            else {
                mailService.sendMail(sentTo,
                        "OPDB job",
                        "Sorry! you job failed, we will finish it\n\ncontact us with ttttttliu@qq.com");
                mailService.sendMail("yd.liu.scu@gmail.com",
                        "OPDB job",
                        sentTo + "\n" + srr_num + "\n" + method + "\n" + kegg_id + "\n\n" + predictResult);
            }

            return "thanks";
        }
        else{
            model.addAttribute("srr_num",srr_num);
            return "result";
        }
    }

}
