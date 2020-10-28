package com.example.rest;

import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

@RestController
public class RepositoryResource {
    
    @RequestMapping
    @GetMapping
    public ResponseEntity<String> getResponse(){
        return ResponseEntity.ok("Hi");
    }
}
