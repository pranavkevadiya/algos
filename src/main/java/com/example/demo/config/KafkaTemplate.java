package com.example.demo.config;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.kafka.support.SendResult;
import org.springframework.stereotype.Component;
import org.springframework.util.concurrent.ListenableFuture;
import org.springframework.util.concurrent.ListenableFutureCallback;

@Component
public class KafkaTemplate {
    
    private Logger LOG = LogManager.getLogger();
    
    @Autowired
    private org.springframework.kafka.core.KafkaTemplate<String, String> kafkaTemplate;
    
    public void sendMessage(String msg) {
        ListenableFuture<SendResult<String, String>> send 
                = kafkaTemplate.send(KafkaTopicConfig.TOPIC_NAME, msg);
        send.addCallback(new ListenableFutureCallback<SendResult<String, String>>() {
            @Override
            public void onFailure(Throwable throwable) {
                LOG.info("Message failed to send", throwable);
    
            }
    
            @Override
            public void onSuccess(SendResult<String, String> stringStringSendResult) {
                LOG.info("Message sent successfully :" + stringStringSendResult.toString());
        
            }
        });
    }
}
