/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.jlab.service.dc;

/**
 *
 * @author ziegler
 */
public class DCTBEngineAI extends DCTBEngine {
    public DCTBEngineAI() {
        super("DCTBAI");
        super.aiAssist = true;
        super._name = "AI";
    }
}
