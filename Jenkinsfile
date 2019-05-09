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
        sh './configure'
        sh 'make'
      }
    }

    stage('Test') {
      steps {
        sh 'make run'
      }
    }
  }

  post { always { cleanWs() } }
}
