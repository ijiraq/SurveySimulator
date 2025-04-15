
REPO = images.canfar.net
PROJECT = uvickbos
DEVNAME = ssim
VERSION = regina2

NAME = $(REPO)/$(PROJECT)/$(DEVNAME)

production: dependencies Dockerfile
	docker build --target deploy --build-arg VERSION=$(VERSION) -t $(NAME):$(VERSION) -f Dockerfile .

deploy: production
	docker push $(NAME):$(VERSION)

dev: dependencies Dockerfile
	docker build --target test --build-arg VERSION=$(VERSION) -t $(NAME):$(VERSION) -f Dockerfile .
	echo "docker run --rm -it -p 8888:8888 $(NAME):$(VERSION) bash"

dependencies: 

init:
	mkdir -p build

.PHONY: clean
clean:
	\rm -rf build
