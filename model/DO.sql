-- MySQL Script generated by MySQL Workbench
-- Mon Apr 23 17:03:17 2018
-- Model: New Model    Version: 1.0
-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

-- -----------------------------------------------------
-- Schema DO
-- -----------------------------------------------------
-- Table storing multiple disease ontology ids for Monarch Initiative

-- -----------------------------------------------------
-- Schema DO
--
-- Table storing multiple disease ontology ids for Monarch Initiative
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `DO` DEFAULT CHARACTER SET utf8 ;
USE `DO` ;

-- -----------------------------------------------------
-- Table `DO`.`DO_entryId`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `DO`.`DO_entryId` (
  `DB` VARCHAR(45) NOT NULL COMMENT 'Name original database/ontology',
  `id` VARCHAR(45) NOT NULL COMMENT 'Disease ontology identifier from Disease Ontology',
  `definition` VARCHAR(250) NULL,
  PRIMARY KEY (`DB`, `id`))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `DO`.`DO_crossId`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `DO`.`DO_crossId` (
  `DB1` VARCHAR(45) NOT NULL COMMENT 'Name database for id1',
  `id1` VARCHAR(45) NOT NULL COMMENT 'disease ontology identifier',
  `DB2` VARCHAR(45) NULL COMMENT 'Name database id2',
  `id2` VARCHAR(45) NULL COMMENT 'Crossreference disease ontology id to id1',
  PRIMARY KEY (`DB1`, `id1`),
  INDEX `fk_crossId_entryId_idx` (`DB1` ASC, `id1` ASC),
  CONSTRAINT `fk_crossId_entryId`
    FOREIGN KEY (`DB1` , `id1`)
    REFERENCES `DO`.`DO_entryId` (`DB` , `id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `DO`.`DO_parentId`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `DO`.`DO_parentId` (
  `DB` VARCHAR(45) NOT NULL COMMENT 'Database for id',
  `id` VARCHAR(45) NOT NULL COMMENT 'Disease ontology identifier from Disease Ontology ',
  `pDB` VARCHAR(45) NULL COMMENT 'Name database for parent id',
  `parent` VARCHAR(45) NULL COMMENT 'Parent ontology for id in Monarch Disease Ontology ',
  PRIMARY KEY (`DB`, `id`),
  INDEX `DB, id_idx` (`pDB` ASC, `parent` ASC),
  CONSTRAINT `fk_table1_entryId1`
    FOREIGN KEY (`DB` , `id`)
    REFERENCES `DO`.`DO_entryId` (`DB` , `id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `DB, id`
    FOREIGN KEY (`pDB` , `parent`)
    REFERENCES `DO`.`DO_entryId` (`DB` , `id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `DO`.`DO_idNames`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `DO`.`DO_idNames` (
  `DB` VARCHAR(45) NOT NULL COMMENT 'Name original database',
  `id` VARCHAR(45) NOT NULL COMMENT 'Disease ontology identifier from Disease Ontology',
  `name` VARCHAR(45) NULL COMMENT 'Term (synonym or label) to describe the disease',
  `canonical` TINYINT NULL COMMENT 'Current label for the entry',
  PRIMARY KEY (`DB`, `id`),
  CONSTRAINT `fk_table1_entryId2`
    FOREIGN KEY (`DB` , `id`)
    REFERENCES `DO`.`DO_entryId` (`DB` , `id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `DO`.`DO_sourceFiles`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `DO`.`DO_sourceFiles` (
  `url` VARCHAR(45) NOT NULL COMMENT 'URL location of source files',
  `current` VARCHAR(45) NULL COMMENT 'Date of the current version of the source files',
  PRIMARY KEY (`url`))
ENGINE = InnoDB;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;