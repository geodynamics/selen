#!groovy

pipeline {
  agent {
    docker {
      image 'geodynamics/selen-buildenv-bionic:latest'
      alwaysPull true
    }
  }

  options {
    timeout(time: 2, unit: 'HOURS')
  }

  stages {
    stage('Build') {
      steps {
        sh 'cd src ; make'
      }
    }

    stage('Test') {
      steps {
        sh 'cp CONFIGS/config.sle.TEST .'
        sh 'cd DATA ; gunzip *R30*.gz'
        sh './sha.exe 30 32 DATA/px-R30.dat DATA/px-lat-R30.dat DATA/sh-R30L32.bin'
        sh 'sh ./make_sle.sh TEST'
      }
    }
  }

  post { always { cleanWs() } }
}
