package edu.scu.service;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.mail.SimpleMailMessage;
import org.springframework.mail.javamail.JavaMailSender;
import org.springframework.stereotype.Service;

@Service
public class MailServiceImpl implements MailService{
    private final JavaMailSender mailSender;

    @Autowired
    public MailServiceImpl(JavaMailSender mailSender) {
        this.mailSender = mailSender;
    }

    @Override
    public String sendMail(String sentTo, String subject, String content){
        String sentFrom = "ttttttliu@qq.com";

        SimpleMailMessage email = new SimpleMailMessage();
        email.setFrom(sentFrom);
        email.setTo(sentTo);
        email.setSubject(subject);
        email.setText(content);

        mailSender.send(email);
        return "OK";
    }
}
