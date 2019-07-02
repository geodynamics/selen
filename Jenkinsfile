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
        sh 'cd DATA ; gunzip *R44*.gz'
        sh './sha.exe 44 128 DATA/px-R44.dat DATA/px-lat-R44.dat DATA/sh-R44L128.bin'
        sh 'sh ./make_sle.sh I6G-R44-L128-I33'
      }
    }
  }

  post { always { cleanWs() } }
}
