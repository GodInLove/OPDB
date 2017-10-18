package edu.scu.service;

public interface MailService {
    public String sendMail(String sentTo, String subject, String content);
}
