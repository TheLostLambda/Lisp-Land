;;Stage 1: Implement Characteristics of Life (2 weeks) : DONE
;;Stage 1.5: Add complementary functions (2 weeks) : DONE
;;Stage 2: Add in realism and research (6+ weeks) : IN PROGRESS
;;Stage 3: Beautify and refine code (TBD) : TODO

;;;; Current Task: Rethink 'Cell-Loc'/add 'Cell-Sft' to avoid stacking

(ql:quickload 'lispbuilder-sdl)
(ql:quickload 'lispbuilder-sdl-gfx)

(defparameter *width* 100)
(defparameter *height* 100)
(defparameter *area* (* *width* *height*))
(defparameter *pixsize* 5) ;;Note + TODO: One "pixel" in the sim is equal to 20um, this is an estimate, revise in stage two.
(defparameter *umperpix* 20) ;;TODO: Use this in stage two to create a correctly scaled measure of size.
(defparameter *cellHP* 150)
(defparameter *cells* nil)
(defparameter *world* nil)
(defparameter *generation* 0)
(defparameter *datafile* "Celldata.db")


(defun btwn (x min max)
  (and (>= x min) (<= x max)))

(defun random-range (min max)
  (+ min (random (- max min) (make-random-state t))))
  
(defun fout (in)
  (format *standard-output* "~a" in)
  (force-output *standard-output*))

(defun remove-nth (n lst) ;;TODO: Make iterative in stage three.
  (if (or (zerop n) (null lst))
      (cdr lst)
      (cons (car lst) (remove-nth (1- n) (cdr lst)))))

(defmacro dotimes-dec (var-and-limit &rest body)
  (let ((var (first var-and-limit))
        (limit (second var-and-limit)))
  `(do ((,var (1- ,limit) (1- ,var)))
       ((<= ,var -1))
    ,@body)))

(defmacro fetch-value (accessor index lst)
  `(getf (nth ,index ,lst) ,accessor))  

(defun seek-propval (accessor value lst)
  (let ((prop-cells '()))
  (mapcar #'(lambda (x)
    (if (equal (getf x accessor) value)
      (append prop-cells x)
      (continue)))
    lst) prop-cells))
    
(defun surround-space (cen)
  (let ((x (car cen)) (y (cdr cen)) (openspcs '()))
  
    (if (seek-propval :POS (cons x (1+ y)) *cells*)
      (continue)
      (append openspcs (cons x (1+ y))))
      
    (if (seek-propval :POS (cons (1+ x) (1+ y)) *cells*)
      (continue)
      (append openspcs (cons (1+ x) (1+ y))))
      
    (if (seek-propval :POS (cons (1+ x) y) *cells*)
      (continue)
      (append openspcs (cons (1+ x) y)))
      
    (if (seek-propval :POS (cons (1+ x) (1- y)) *cells*)
      (continue)
      (append openspcs (cons (1+ x) (1- y))))
      
    (if (seek-propval :POS (cons x (1- y)) *cells*)
      (continue)
      (append openspcs (cons x (1- y))))
      
    (if (seek-propval :POS (cons (1- x) (1- y)) *cells*)
      (continue)
      (append openspcs (cons (1- x) (1- y))))
      
    (if (seek-propval :POS (cons (1- x) y) *cells*)
      (continue)
      (append openspcs (cons (1- x) y)))
      
    (if (seek-propval :POS (cons (1- x) (1+ y)) *cells*)
      (continue)
      (append openspcs (cons (1- x) (1+ y))))
      openspcs))
  
(defun fetch-props (lst)
  (let ((nlst nil))
    (do ((i 0 (+ i 2)))
        ((>= i (length lst))) ;;TODO: Modify to deal with non-sequences in stage 3.
      (push (nth i lst) nlst)) (reverse nlst)))  

(defun new-cell (POS ATP NA AA FA G O2 CO2 DNA)
  (push (list :POS POS :ATP ATP :NA NA 
         :AA AA :FA FA :G G :O2 O2 :CO2 CO2 :DNA DNA :HP *cellHP*) *cells*))

(defun new-world (TEMP PH NAC AAC FAC GC O2C CO2C RAD BPV LUX)
  (setf *world* (list :TEMP TEMP :PH PH :NAC NAC :AAC AAC :FAC FAC :GC GC
                      :O2C O2C :CO2C CO2C :RAD RAD :BPV BPV :LUX LUX)))
                
(defun parse-gene (gene)
  (let ((gene-val 0))
    (dotimes (i (length gene))
      (setf gene-val (+ (* (expt 2 i) (nth i (reverse gene))) gene-val))) gene-val))
      
(defun parse-DNA (DNA-seq)
  (let ((DNA-val nil))
    (dotimes (i (length DNA-seq))
      (push (parse-gene (nth i DNA-seq)) DNA-val)) (reverse DNA-val)))

(defun dump-ci ()
  (dolist (entry *cells*)
    (format t "~%~{~a:~10t~a~%~}~%" entry)))

(defun save-sim (filename)
  (with-open-file (out filename
                   :direction :output
                   :if-exists :supersede)
    (with-standard-io-syntax
      (print *cells* out)
      (print *world* out))))

(defun load-sim (filename)
  (with-open-file (in filename)
    (with-standard-io-syntax
      (setf *cells* (read in))
      (setf *world* (read in)))))
      
(defun rand-mutate (DNA-seq)
  (let ((CDNA-seq (copy-tree DNA-seq)))
    (let ((gene (random-range 0 (length CDNA-seq))))
      (let ((base (random-range 0 (length (nth gene CDNA-seq)))))
        (cond ((= (nth base (nth gene CDNA-seq)) 0) (setf (nth base (nth gene CDNA-seq)) 1))
              (t (setf (nth base (nth gene CDNA-seq)) 0))) CDNA-seq))))
            
(defun mutate-times (DNA-seq times)
  (dotimes (i times)
    (rand-mutate DNA-seq)))
  
(defun close-sim ()
  (format t "Saving simulation state and exiting...~%")
  (save-sim *datafile*)
  (exit))
  
(defun init-sim ()
  (format t "Loading simulation state and starting...~%")
  (load "GeneR.lisp")
  (when (probe-file *datafile*)
    (load-sim *datafile*))
  (format t "Done."))
  
(defun Cell-Apo (celli) ;;Note: Function for cellular death and anything HP related.
  (if (<= (fetch-value :HP celli *cells*) 0)
      (progn ;;TODO: Once 'Cell-Env' can expel molecules, replace this code with a 'Cell-Env' call.
      (incf (getf *world* :NAC) (fetch-value :NA celli *cells*))
      (incf (getf *world* :AAC) (fetch-value :AA celli *cells*))
      (incf (getf *world* :FAC) (fetch-value :FA celli *cells*))
      (incf (getf *world* :GC) (fetch-value :G celli *cells*))
      (incf (getf *world* :O2C) (fetch-value :O2 celli *cells*))
      (incf (getf *world* :CO2C) (fetch-value :CO2 celli *cells*))
      (setf *cells* (remove-nth celli *cells*)))
      ;;TODO: Put other HP effecting conditions here...
      (decf (fetch-value :HP celli *cells*) 1) ;;Note: Cell aging condition.
      ))  

(defun Cell-Vir (celli) ;;Note + TODO: Viruses are pretty crazy, there is a lot to do here, revise in stage two...
  (let ((BPVP nil) (Chance nil) (CPS nil) (NBPVP nil)) 
    (setf BPVP (* (/ (getf *world* :BPV) *area*) 100))
    (setf Chance (random-range 1 101))
    (setf CPS (truncate (/ BPVP 100)))
    (setf NBPVP (- BPVP (* CPS 100)))
    (when (<= NBPVP Chance) (incf CPS 1))
    
    (dotimes (i CPS) 
      (when (= (random-range 1 (+ (GeneR celli 'perm) 1)) 1)
        (decf (fetch-value :HP celli *cells*) 10)
        (incf (getf *world* :BPV) 1)))))
          
(defun Cell-Env (celli) ;;Note + TODO: This is scientifically flawed, revise in stage two...
  (let ((prop nil) (cprop nil) (propcont nil) (cpropcont nil))
    (do ((i 0 (1+ i)))
        ((>= i 6)) ;;Note + TODO: 6 is for the four main macromolecules plus CO2 and O2, make more portable later...
      
      (setf prop (nth i '(:NAC :AAC :FAC :GC :O2C :CO2C))) ;;TODO: Make more portable later...
      (setf cprop (nth i '(:NA :AA :FA :G :O2 :CO2))) ;;TODO: Make more portable later...
      (setf propcont (getf *world* prop))
      (setf cpropcont (fetch-value cprop celli *cells*))
      
      ;;NOTE + TODO: This is a hiding place for bugs, review this.
      (cond ((> cpropcont (/ propcont *area*))
              (let ((CPSD (truncate (- cpropcont (/ propcont *area*)))))
                (dotimes (i CPSD)
                  (when (= (random-range 1 (+ (GeneR celli 'perm) 1)) 1)
                    (decf (fetch-value cprop celli *cells*) 1)
                    (incf (getf *world* prop))))))
            ((< cpropcont (/ propcont *area*))
              (let ((PropVal nil) (Chance nil) (CPS nil) (NPropVal nil)) 
				(setf PropVal (* (/ (getf *world* prop) *area*) 100))
				(setf Chance (random-range 1 101))
				(setf CPS (truncate (/ PropVal 100)))
				(setf NPropVal (- PropVal (* CPS 100)))
				(when (<= NPropVal Chance) (incf CPS 1))
    
				(dotimes (i CPS) 
				  (when (= (random-range 1 (+ (GeneR celli 'perm) 1)) 1)
					(decf (getf *world* prop) 1)
					(incf (fetch-value cprop celli *cells*) 1)))))
            (t (continue))))))
            
(defun Cell-Pho (celli) ;;TODO: Incorperate Lux into the photosynthesis process.
  (when (> (GeneR celli 'chlorop) 1)
    (cond ((> 1 0)
      (dotimes (i 4) ;;Dummy Value: 4 is a place-holder for the chloroplast count.
        (when (and (> (fetch-value :CO2 celli *cells*) 6) (> (fetch-value :O2 celli *cells*) 6))
        (decf (fetch-value :CO2 celli *cells*) 6) (incf (fetch-value :O2 celli *cells*) 6)
        (incf (fetch-value :G celli *cells*) 1)))) ;;Note: Photosynthesis 6 CO2 + 6 H2O → C6H12O6 + 6 O2. For now, Water is ignored
      (t (format t "If you see this message, the laws of science have broken down.~%Now is the time for panic")))))
  
(defun Cell-Met (celli) ;;TODO: Add a better way to regulate the ATP synthesis
  (when (< (fetch-value :ATP celli *cells*) 1500) ;;Dummy Value: 1500 is the place-holder for max APT value.
    (cond ((> 1 0) ;;Note + TODO: This is the condition for mitochondrial respiration, revise during stage two...
      (dotimes (i 2) ;;Dummy Value: 2 is a place-holder for the mitochondron count.
       (cond ((and (>= (fetch-value :O2 celli *cells*) 6) (>= (fetch-value :G celli *cells*) 1)) ;;Note: Aerobic Respiration C6H12O6 + 6 O2 → 6 CO2 + 6 H2O + 38 ATP. For now, Water is ignored
         (decf (fetch-value :G celli *cells*) 1) (decf (fetch-value :O2 celli *cells*) 6)
         (incf (fetch-value :CO2 celli *cells*) 6) (incf (fetch-value :ATP celli *cells*) 38))
           
             ((and (< (fetch-value :O2 celli *cells*) 6) (>= (fetch-value :G celli *cells*) 1)) ;;Note: Anaerobic Respiration C6H12O6 → 2 CO2 + 2 C2H5OH + 2 ATP. For now Ethanol is ignored
         (decf (fetch-value :G celli *cells*) 1) (incf (fetch-value :CO2 celli *cells*) 2)
         (incf (fetch-value :ATP celli *cells*) 2)))))
      ;;TODO: This is the condition for ATP synthesis via a proton gradient.   
      (t (format t "If you see this message, the laws of science have broken down.~%Now is the time for panic")))))  
  
;;TODO: Prevent cells from overlapping  
(defun Cell-Loc (celli &optional (dir (random 33 (make-random-state t)))) ;;Note + Dummy Value: This is the function for cellular locomotion. 
  (let ((x (car (fetch-value :POS celli *cells*))) (y (cdr (fetch-value :POS celli *cells*)))) ;;TODO: Find a more realistic and controlled method of locomotion.
    (when (and (>= (fetch-value :ATP celli *cells*) 25) (<= dir 8)) ;;Dummy Value: 25 is ATP cost for movement.
	  (decf (fetch-value :ATP celli *cells*) 25) ;;Dummy Value: 25 is ATP cost for movement.
	  (cond ((= dir 1) (setf (fetch-value :POS celli *cells*) (cons x (mod (1+ y) *height*))))
	        ((= dir 2) (setf (fetch-value :POS celli *cells*) (cons (mod (1+ x) *width*) y))
		               (setf (fetch-value :POS celli *cells*) (cons x (mod (1+ y) *height*))))
		    ((= dir 3) (setf (fetch-value :POS celli *cells*) (cons (mod (1+ x) *width*) y)))
		    ((= dir 4) (setf (fetch-value :POS celli *cells*) (cons (mod (1+ x) *width*) y))
		               (setf (fetch-value :POS celli *cells*) (cons x (mod (1- y) *height*))))
		    ((= dir 5) (setf (fetch-value :POS celli *cells*) (cons x (mod (1- y) *height*))))
		    ((= dir 6) (setf (fetch-value :POS celli *cells*) (cons (mod (1- x) *width*) y))
		               (setf (fetch-value :POS celli *cells*) (cons x (mod (1- y) *height*))))
		    ((= dir 7) (setf (fetch-value :POS celli *cells*) (cons (mod (1- x) *width*) y)))
		    ((= dir 8) (setf (fetch-value :POS celli *cells*) (cons (mod (1- x) *width*) y))
		               (setf (fetch-value :POS celli *cells*) (cons x (mod (1+ y) *height*))))))))
  
(defun Cell-Mut (celli)
  (when (<= (random-range 1 101) (getf *world* :RAD))
    (setf (fetch-value :DNA celli *cells*) (rand-mutate (fetch-value :DNA celli *cells*)))))

;;NOTE + TODO: Truncate here is breaking the conservation law. Some molicules are just destroyed. Fix this.  
(defun Cell-Rep (celli)
  (when (>= (fetch-value :ATP celli *cells*) 1500)
    (let ((POS (fetch-value :POS celli *cells*))
          (ATP (truncate (/ (- (fetch-value :ATP celli *cells*) 1000) 2)))
          (NA (truncate (/ (fetch-value :NA celli *cells*) 2)))
          (AA (truncate (/ (fetch-value :AA celli *cells*) 2)))
          (FA (truncate (/ (fetch-value :FA celli *cells*) 2)))
          (G (truncate (/ (fetch-value :G celli *cells*) 2)))
          (O2 (truncate (/ (fetch-value :O2 celli *cells*) 2)))
          (CO2 (truncate (/ (fetch-value :CO2 celli *cells*) 2)))
          (DNA (fetch-value :DNA celli *cells*)))
  
      (setf (fetch-value :ATP celli *cells*) ATP)
      (setf (fetch-value :NA celli *cells*) NA)
      (setf (fetch-value :AA celli *cells*) AA)
      (setf (fetch-value :FA celli *cells*) FA)
      (setf (fetch-value :G celli *cells*) G)
      (setf (fetch-value :O2 celli *cells*) O2)
      (setf (fetch-value :CO2 celli *cells*) CO2)
      (new-cell POS ATP NA AA FA G O2 CO2 DNA)
      (Cell-Loc 0 (random 9 (make-random-state t))))))
      
(defun Cell-Sft (celli) ;;Bug: This hangs the program
  (let ((cellpos (fetch-value :POS celli *cells*)))
    (if (= (length (seek-propval :POS cellpos *cells*)) 1)
      (continue)
      (let ((sur-spc (surround-space cellpos)))
        (if (equal sur-spc nil)
          (progn (Cell-Loc celli (random 9 (make-random-state t))) (Cell-Sft celli))
          (setf (fetch-value :POS celli *cells*) (nth (random-range 0 (length sur-spc)) sur-spc)))))))
  
(defun next-tick ()
  (incf *generation* 1)

  (dotimes-dec (i (length *cells*))
    (Cell-Env i))
  (dotimes-dec (i (length *cells*))
	(Cell-Pho i))
  (dotimes-dec (i (length *cells*))
    (Cell-Met i))
  (dotimes-dec (i (length *cells*))
	(Cell-Loc i))
  (dotimes-dec (i (length *cells*))
	(Cell-Mut i))
  (dotimes-dec (i (length *cells*))
	(Cell-Vir i))
  (dotimes-dec (i (length *cells*))
	(Cell-Rep i))
  (dotimes-dec (i (length *cells*))
	(Cell-Apo i))
  (dotimes-dec (i (length *cells*))
	(Cell-Sft i))
  
  (fout *generation*)
  (fresh-line))

(defun autoplay-sim (&optional (delay 0))
  (sdl:with-init ()
  (sdl:window (* *pixsize* *width*) (* *pixsize* *height*) :title-caption "Lisp Land")
    (setf (sdl:frame-rate) 60)
    (sdl:with-events ()
      (:quit-event () t)
      (:idle ()
        (sleep delay)
        (next-tick)
        (sdl:clear-display sdl:*black*)
        (dotimes (celli (length *cells*))
          (let ((POS (fetch-value :POS celli *cells*)))
                  (sdl-gfx:draw-box (sdl:rectangle :x (* (car POS) *pixsize*) :y (* (cdr POS) *pixsize*)
                                    :w *pixsize* :h *pixsize*) :color (if (> (GeneR celli 'chlorop) 1) sdl:*green* sdl:*white*))))
        (sdl:update-display)))))

(init-sim)
